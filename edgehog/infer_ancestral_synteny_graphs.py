#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import itertools
import networkx as nx
import numpy as np
import pyham


def print_graph_info(genome, graph):
    average_degree = np.mean([i[1] for i in graph.degree()])
    print('%s graph: nodes=%d; edges=%d; connected_components=%d; avg_degree=%.2f' % 
                (genome.name, graph.number_of_nodes(), graph.number_of_edges(), 
                 nx.number_connected_components(graph), average_degree))


####################################################################################
# GOAL:
#   the transitive graph is the object that allows to propagate edges from a child 
#   to its parent node on the phylogeny

# WHAT IT DOES:
#   Modifies the graph such that genes that do not have an ancestral gene
#   (rootHOG or gene without orthologs) are removed and their adjacent nodes are reconnected
#   (as long as number of adjacent genes without parent <= max_gap_length)

# rootHOG OF DEGREE > 2:
#   if max_rootHog_neighbors_for_reconnexion == 2:
#     the reconnexion is allowed only if the rootHOG has stricly two neighbors: 
#       A-b-C --> A-C while A-b-C --> A; B; C (break point: A, B and C disconnected)
#                             |
#                             D
#   if max_rootHog_neighbors_for_reconnexion > 2:
#     all unique pairs of neighbors of a rootHOG will be reconnected. All neighbors are then forming a clique.
#     A-b-C --> A-C 
#       |       \ /
#       D        D
#   the top-down phase is assumed to do the job for deciphering which edge of the clique is actually true

# MAX GAP LENGTH:
#   This parameter can be seen as the theoritical maximal number of consecutive novel genes 
#   that can emerge / be inserted between two older genes. For instance:
#     if max_gap_length == 2: 
#        the probabilistic A-b-c-D-E-f-g-h-I-J graph is turn into A-D-E ; I-J (2 CCs)
#     whereas if max_gap_length == 3:
#        the probabilistic A-b-c-D-E-f-g-h-I-J is turn into A-D-E-I-J (1 CC)

# NOTE REGARDING MAX GAP LENGTH FOR PROKARYOTES:
#   in prokaryotes where insertion of exogenic genomic fragments is frequent,
#   (e.g. through foreign DNA uptake via competence, integrative conjugative element (ICE) insertion, prophage insertion etc)
#   it might be relevant to increase the max_gap_length
#

# WEIGHT ASSIGNMENT:
#    upper case: gene present in ancestor / lower case: gene absent in ancestor
#    e.g      A-C    A-C-E  A-E
#    weight =  6      1 1    1
#               \   /       /
#               A-c-E      /
#    weight =    7 1      /
#                  \     /
#                    A-E 
#    weight =         8
####################################################################################     
def build_transitive_graph(graph, max_gap_length = 2):
    real_edge_to_transitive_edge = dict()
    reconnected_edge_to_real_edges = dict()
    cnt_removed_nodes = 0
    max_rootHog_neighbors_for_reconnexion = 2
    gap_length = dict()
    genes = list(graph.nodes)
    for gene in genes:
        if gene.parent is None:
            neighbors = list(graph.adj[gene]) 
            if len(neighbors) > 1 and len(neighbors) <= max_rootHog_neighbors_for_reconnexion:
                # create a clique between neighbors (the top down phase will naturally prune unsupported edges in the clique)
                for left_gene, right_gene in itertools.combinations(neighbors, 2):
                    gap = 1
                    # increase gap if left_gene and current_gene have already been reconnected
                    if (left_gene, gene) in gap_length:
                        gap += gap_length[(left_gene, gene)]
                    if (gene, right_gene) in gap_length:
                        gap += gap_length[(gene, right_gene)]
                    if gap <= max_gap_length:
                        if not graph.has_edge(left_gene, right_gene):
                             graph.add_edge(left_gene, right_gene, 
                                            weight = min(graph[left_gene][gene]["weight"], graph[gene][right_gene]["weight"]),
                                            children = list(set(graph[left_gene][gene]["children"] + graph[gene][right_gene]["children"])),
                                            extant_descendants = list(set(graph[left_gene][gene]["extant_descendants"] + graph[gene][right_gene]["extant_descendants"])))
                             # if gene that will be removed has already been reconnected:
                             if (left_gene, gene) in reconnected_edge_to_real_edges:
                                 for e in reconnected_edge_to_real_edges[(left_gene, gene)]:
                                     real_edge_to_transitive_edge[e] = (left_gene.parent, right_gene.parent)
                                     real_edge_to_transitive_edge[e[::-1]] = (left_gene.parent, right_gene.parent)
                             if (gene, right_gene) in reconnected_edge_to_real_edges:
                                 for e in reconnected_edge_to_real_edges[(gene, right_gene)]:
                                     real_edge_to_transitive_edge[e] = (left_gene.parent, right_gene.parent)
                                     real_edge_to_transitive_edge[e[::-1]] = (left_gene.parent, right_gene.parent)
                             real_edge_to_transitive_edge[(left_gene, gene)] = (left_gene.parent, right_gene.parent)
                             real_edge_to_transitive_edge[(gene, left_gene)] = (left_gene.parent, right_gene.parent)
                             real_edge_to_transitive_edge[(right_gene, gene)] = (left_gene.parent, right_gene.parent)
                             real_edge_to_transitive_edge[(gene, right_gene)] = (left_gene.parent, right_gene.parent)
                             # the following temp dictionary will be used to update the conversion real edge to transitive edge 
                             # if a previously assigned transitive edge ends up being trimmed as well (consecutive removals (gap > 1))
                             reconnected_edge_to_real_edges[(left_gene, right_gene)] = [(left_gene, gene), (gene, right_gene)]
                             reconnected_edge_to_real_edges[(right_gene, left_gene)] = [(left_gene, gene), (gene, right_gene)]
                        # when reconnecting two nodes, the gap_length is updated
                        gap_length[(left_gene, right_gene)] = gap_length[(right_gene, left_gene)] = gap
                graph.remove_node(gene)
                cnt_removed_nodes += 1
    # print("      - removed %d nodes that have no parent hog" % cnt_removed_nodes)
    return real_edge_to_transitive_edge

   
####################################################################################         
# BOTTOM-UP PHASE

# GOAL:
#   Propagate count of adjacencies observed between two genes 
#   from children nodes to their parent node on the species tree 
#   (as long as the parent node has the two ancestral genes)

# RATIONALE:
#   by parsimony, the more extant genes belonging to HOGa and HOGb 
#   are found adjacent to each other in extant genomes, 
#   the more likely ancestral genes belonging to HOGa and HOGb 
#   were adjacent to each other in the ancestor of these genomes
#################################################################################### 
def leaves_to_root_synteny_propagation(ham, max_gaps):
    print('###################################')
    print('Bottom-up_phase: propagate count of adjacencies observed between two genes from children nodes to their parent node on the species tree '
          '(as long as the parent node has the two ancestral genes)')
    nb_ancestors = ham.taxonomy.tree.n_internal_nodes
    anc_counter = 0
    for tree_node in ham.taxonomy.tree.traverse('postorder'):
        try:
            genome = tree_node.genome
        except AttributeError:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if isinstance(genome, pyham.AncestralGenome):
            anc_counter += 1
            print('processing ancestral genome %d/%d: %s'  % (anc_counter, nb_ancestors, genome.name))
            graph = nx.Graph()
            graph.add_nodes_from(genome.genes)
                
            for child in tree_node.children:
                # add a node in the parent graph, whenever a newly visited child gene with a parent is encountered
                # this is at the basis of ancestral gene content reconstruction and allows for singletons in ancestral graphs
                child_working_graph = child.bottom_up_synteny.copy()
                for gene in child_working_graph:
                    if gene.parent is not None and not graph.has_node(gene.parent):
                        graph.add_node(gene.parent)
                        
                
                # For the child of current ancestor, a transitivity graph allowing comparison with ancestral gene content is created:
                # genes that do not have a parent hog are removed and their neighbors are connected
                # in a working copy that is used
                # to propagate the synteny graph to the higher level.
                #
                # The trimmed working copy is eventually stored as a transitive graph and will be used subsequently in the top-down phase
                real_edge_to_transitive_edge = build_transitive_graph(child_working_graph, max_gaps)
                child_transitive_graph = nx.Graph()
                for gene, adjacent_gene, weight in child_working_graph.edges.data('weight', default=1):
                    if(gene.parent is not None and adjacent_gene.parent is not None and gene.parent != adjacent_gene.parent):
                        assert gene.parent in graph.nodes
                        assert adjacent_gene.parent in graph.nodes
                        ancestral_edge = (gene.parent, adjacent_gene.parent)
                        
                        if isinstance(child.genome, pyham.ExtantGenome):
                            extant_descendants = [child]
                        else:
                            extant_descendants = child_working_graph[gene][adjacent_gene]['extant_descendants']
                            
                        # update edge weight if edge already exists in the graph (e.g. when second children is visited)
                        if graph.has_edge(*ancestral_edge):
                            graph[ancestral_edge[0]][ancestral_edge[1]]['weight'] += weight
                            graph[ancestral_edge[0]][ancestral_edge[1]]['children'].append(child)
                            graph[ancestral_edge[0]][ancestral_edge[1]]['extant_descendants'] += extant_descendants
                        else:
                            # create edge otherwise, with the same weight as in current supporting child
                            graph.add_edge(*ancestral_edge, weight = weight, children = [child], extant_descendants = extant_descendants)
                                
                        if not gene.parent in child_transitive_graph.nodes:
                            child_transitive_graph.add_node(gene.parent, real_gene = gene)
                        if not adjacent_gene.parent in child_transitive_graph.nodes:
                            child_transitive_graph.add_node(adjacent_gene.parent, real_gene = adjacent_gene)
                        if child_transitive_graph.has_edge(*ancestral_edge):
                            child_transitive_graph[gene.parent][adjacent_gene.parent]['weight'] += weight
                        else:
                            child_transitive_graph.add_edge(*ancestral_edge, weight = weight)
                        
                child.add_feature('transitive_synteny', child_transitive_graph)
                child.add_feature('real_edge_to_transitive_edge', real_edge_to_transitive_edge)
            
            tree_node.add_feature('bottom_up_synteny', graph)
            
####################################################################################         
# TOP-DOWN PHASE
 
# GOAL:
#    1) Infer in which node of the tree an adjacency between two genes has emerged 
#       (The LCA is identified as the last internal node for which the adjacency is observed in two direct children)
#    2) Prune this adjacency in all ancestors predating this node 
#       (might have been propagated in the past if the two genes happened to be present in the ancestors of this node) 

# RATIONALE:
#   an adjacency between two genes is unlikely in an ancestral genome 
#   if this adjacency is observed only in one child and not in the outgroup
#   (in this case, the LCA for this adjacency should be in the child lineage)    
####################################################################################           
def root_to_leaves_edge_trimming(ham):
    print('###################################')
    print('Top-down phase: prune any adjacency propagated before the last ancestor in which this adjacency is inferred to have emerged ')
    nb_ancestors = ham.taxonomy.tree.n_internal_nodes
    anc_counter = 0
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        try:
            genome = tree_node.genome
        except AttributeError:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if isinstance(genome, pyham.AncestralGenome):
            anc_counter += 1
            print('processing ancestral genome %d/%d: %s'  % (anc_counter, nb_ancestors, genome.name))
            edges_to_remove = list()
            graph = tree_node.bottom_up_synteny.copy()
            k = 0
            for gene, adjacent_gene, weight in graph.edges.data('weight', default=1):
                k += 1
                remove_edge = True
                for child in tree_node.children:
                    child_transitive_graph = child.transitive_synteny
                    try:
                        child_weight = child_transitive_graph[gene][adjacent_gene]['weight']
                    except:
                        child_weight = 0
                    # if at least two children contribute to the weight of the edge,
                    # that means the edge was likely present in the current node
                    if child_weight > 0 and child_weight < weight:
                        remove_edge = False
                        graph[gene][adjacent_gene]['nb_internal_nodes_from_ancestor_with_updated_weight'] = 0
                        break
                # if edge is only present in one children lineage, check if it is supported by outgroup (retain) or not (discard)
                if remove_edge and tree_node.up:
                    try:
                        # either there was some reconnection between genes because of insertion in descendant
                        transitive_edge = tree_node.real_edge_to_transitive_edge[(gene, adjacent_gene)]
                    except:
                        # or the real edge in current genome was also present in ancestor
                        transitive_edge = (gene.parent, adjacent_gene.parent)
                    ancestral_graph = tree_node.up.top_down_synteny
                    graph[gene][adjacent_gene]['nb_internal_nodes_from_ancestor_with_updated_weight'] = 1
                    if transitive_edge[0] == transitive_edge[1] and transitive_edge[0] in ancestral_graph.nodes:
                        remove_edge = False
                        # when adjacency between paralogs has been merged in the internal node just before the duplication
                        # then the edge A1-A2 in descendant is simply justified by the presence of A in the ancestor 
                    elif ancestral_graph.has_edge(*transitive_edge):
                        ancestral_weight = ancestral_graph[transitive_edge[0]][transitive_edge[1]]['weight']
                        remove_edge = False
                        if ancestral_weight <= weight:
                            # that means that the edge is not present in the direct outgroup but is present
                            # in extant species connected to an even more ancestral node
                            # in this old ancestral node, the weight of the edge > edge in current node 
                            # we count the number of internal nodes separating this old ancestor from current node
                            # if counter high: possible HGT or edge lost once and reemerged afterwards (convergent evolution)
                            if 'nb_internal_nodes_from_ancestor_with_updated_weight' in ancestral_graph[transitive_edge[0]][transitive_edge[1]]:
                                graph[gene][adjacent_gene]['nb_internal_nodes_from_ancestor_with_updated_weight'] += ancestral_graph[transitive_edge[0]][transitive_edge[1]]['nb_internal_nodes_from_ancestor_with_updated_weight']
                if remove_edge:
                    edges_to_remove.append((gene, adjacent_gene))
            graph.remove_edges_from(edges_to_remove)
            tree_node.add_feature('top_down_synteny', graph)
            

def compute_neighbor_cost(graph, node, previous_node, max_path_len = 500):
    cost = graph[previous_node][node]['weight']
    path_len = 1
    while graph.degree(node) == 2 and path_len < max_path_len:
        next_node = [neighbor for neighbor in graph.neighbors(node) if neighbor != previous_node][0]
        previous_node = node
        node = next_node
        cost += graph[previous_node][node]['weight']
        path_len += 1
    return cost

####################################################################################
# LINEARIZATION

# GOAL:
#   for each ancestral synteny graph,
#   transform each connected component with conflicting nodes (degree>2) into the most likely linear path(s)  

# RATIONALE:
#   as a real genome is linear/circular, a gene cannot have more than two direct neighbors

# WHAT IT DOES:
#   Resolve conflicting nodes by pruning their edges to worse neighbors until degree becomes 2.
#   Conflicting nodes are resolved by ascending number of degree.
#   In case of weight equality of between several neighbors, 
#   the neighbor that minimizes the weight of the resulting connected component is chosen first for pruning

# STEPS:
#   step 1: sort nodes by ascending number of degree 
#     for each degree > 2:
#       sort corresponding nodes by descending number of the sum of the weights to their direct neighbors
#   step 2: for each conflicting node in the sorted list:
#       for each neighbor: 
#         get the path from the source node that passes by the current neighbor and stops at the first node encountered with degree != 2
#       for each neighbor path:
#         get the sum of weights
#       sort neighbors by ascending number of this sum 
#       (useful if subsequent equality of minimal weight to multiple neighbors (we want to maximize the size of the CC))
#   step 3: resort neighbors by ascending number of the weight of the edge to current_node
#       cut edge towards neighbor with minimal weight until degree of current node == 2
#       if same minimal weight, neighbor with the minimal sum of weights in the path that ends with a node of degree != 2 will be favored due to prior sorting
####################################################################################
def linearize_graphs(ham):
    print('###################################')
    print('Linearization: infer the most likely ancestral genome ...')

    nb_ancestors = ham.taxonomy.tree.n_internal_nodes
    anc_counter = 0
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        try:
            genome = tree_node.genome
        except AttributeError:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if isinstance(genome, pyham.AncestralGenome):
            anc_counter += 1
            print('processing ancestral genome %d/%d: %s'  % (anc_counter, nb_ancestors, genome.name))
            graph = tree_node.top_down_synteny.copy()
        
            conflicting_nodes = [(n, d) for (n, d) in sorted(graph.degree(), key=lambda pair: pair[1]) if d > 2]
            sum_weights = list()
            for node, degree in conflicting_nodes:
                sum_weights.append(sum([graph[node][neighbor]['weight'] for neighbor in graph.neighbors(node)]))
            node_order = np.argsort(sum_weights)[::-1].tolist()
            conflicting_nodes = [conflicting_nodes[i] for i in node_order]
            
            removed_edges = list()
            removed_edge_attributes = dict()
            
            degrees = sorted(set([d for (n, d) in conflicting_nodes]))
            for degree in degrees:

                nodes = [n for (n, d) in conflicting_nodes if d == degree]
                for node in nodes:
                    if graph.degree(node) > 2:
                        neighbors_cost = list()
                        neighbors = list()
                        for neighbor in graph.neighbors(node):
                            neighbors.append(neighbor)
                        # sort by degree (in case of weight equality, we want as a secondary criterion to chose the neighbor with the highest degree)   
                        neighbors_degree = [graph.degree(neighbor) for neighbor in neighbors]
                        degree_order = np.argsort(neighbors_degree).tolist()
                        neighbors = [neighbors[i] for i in degree_order]
                        # sort by path cost (in case of weight equality, we want as a primary criterion to chose the neighbor that maximizes CC weight)   
                        # neighbors_cost = [compute_neighbor_cost(graph, neighbor, node, 0) for neighbor in neighbors]
                        neighbors_cost = [compute_neighbor_cost(graph, neighbor, node) for neighbor in neighbors]
                        path_order = np.argsort(neighbors_cost).tolist()
                        neighbors = [neighbors[i] for i in path_order]
                        weights = [graph[node][neighbor]['weight'] for neighbor in neighbors]
                        trim_order = np.argsort(weights).tolist()
                        
                        for i in range(0, len(trim_order)-2):
                            node_a, node_b = node, neighbors[trim_order[i]]
                            removed_edges.append((node_a, node_b))
                            removed_edge_attributes[(node_a, node_b)] = graph[node_a][node_b] 
                            graph.remove_edge(node_a, node_b)
                
                # sometimes, the best neighbor for a conflicting node is not necessarily reciprocal, which can result in a few artefactual edge removals
                # the following code is to correct the issue 
                #     A---4---B---4---C                         A---4---B---4---C                    A---4---B---4---C
                #         2 /   \ 2          -Linearization->                       -Correction->
                # W---4---X---2---Y---4---Z                   W---4---X   Y---4---Z              W---4---X---2---Y---4---Z
                # X-Y is flagged as a potential edge to reconnect because degree of X and Y is < 2
                # (in the above case, B can be for instance a transposase that mess things up)
                idx_edges_to_reconnect = [i for i in range(0, len(removed_edges)) if graph.degree(removed_edges[i][0]) < 2 and graph.degree(removed_edges[i][1]) < 2]
                if len(idx_edges_to_reconnect) > 0:
                    edges_to_reconnect = [removed_edges[i] for i in idx_edges_to_reconnect]
                    weights = [removed_edge_attributes[e]["weight"] for e in edges_to_reconnect]
                    costs = list()
                    for e in edges_to_reconnect:
                        cost = 0
                        for node in e:
                            if graph.degree(node) == 1:
                                cost += compute_neighbor_cost(graph, next(graph.neighbors(node)), node)
                        costs.append(cost)
                    cost_order = np.argsort(costs)[::-1].tolist()
                    edges_to_reconnect = [edges_to_reconnect[i] for i in cost_order]
                    weights = [weights[i] for i in cost_order]
                    idx_edges_to_reconnect = [idx_edges_to_reconnect[i] for i in cost_order]
                    weight_order = np.argsort(weights)[::-1].tolist()
                    idx_to_pop = list()
                    for i in range(0, len(weight_order)):
                        edge = edges_to_reconnect[weight_order[i]]
                        if graph.degree(edge[0]) < 2 and graph.degree(edge[1]) < 2 and not graph.has_edge(*edge):
                            graph.add_edge(edge[0], edge[1], **removed_edge_attributes[(edge[0], edge[1])])
                            idx_to_pop.append(idx_edges_to_reconnect[weight_order[i]])
                    idx_to_pop = sorted(idx_to_pop, reverse = True)
                    for idx in idx_to_pop:
                        if idx < len(removed_edges):
                            removed_edges.pop(idx)
                    
            tree_node.add_feature('linear_synteny', graph)
            
            

        
        