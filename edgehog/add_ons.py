#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pyham


def date_edges(ham):
    print('###################################')
    print('Predicting age of extant edges ...')
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        try:
            genome = tree_node.genome
        except AttributeError:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if isinstance(genome, pyham.AncestralGenome):
            graph = tree_node.top_down_synteny
        elif isinstance(genome, pyham.ExtantGenome):
            graph = tree_node.bottom_up_synteny
        if tree_node.up:
            parent_graph = tree_node.up.top_down_synteny
        for gene, adjacent_gene, edge_attributes in graph.edges.data():
            if not tree_node.up:
                graph[gene][adjacent_gene]["lca_node"] = tree_node
            else:
                if gene.parent is None or adjacent_gene.parent is None:
                    graph[gene][adjacent_gene]["lca_node"] = tree_node
                else:
                    parent_edge = (gene.parent, adjacent_gene.parent)
                    if parent_graph.has_edge(*parent_edge):
                        graph[gene][adjacent_gene]["lca_node"] = parent_graph[gene.parent][adjacent_gene.parent]["lca_node"]
                    else:
                        graph[gene][adjacent_gene]["lca_node"] = tree_node
            graph[gene][adjacent_gene]["lca"] = graph[gene][adjacent_gene]["lca_node"].genome.name
            if isinstance(genome, pyham.ExtantGenome):
                graph[gene][adjacent_gene]["age"] = 1 - graph[gene][adjacent_gene]["lca_node"].red


def phylostratify(ham):
    print('###################################')
    print('Phylostratigraphy ...')
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        try:
            genome = tree_node.genome
        except AttributeError:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if tree_node.up:
            if isinstance(genome, pyham.AncestralGenome):
                graph = tree_node.top_down_synteny
            elif isinstance(genome, pyham.ExtantGenome):
                graph = tree_node.bottom_up_synteny
            parent_graph = tree_node.up.top_down_synteny
            
            parent_occurrence = dict()
            for gene in graph:
                if gene.parent is not None:
                    if not gene.parent in parent_occurrence:
                        parent_occurrence[gene.parent] = 0
                    parent_occurrence[gene.parent] += 1
                    
            stratigraphy = dict()
            for key in ["retained_edges", "duplicated_edges", "gained_edges_due_to_gene_gain",
                        "gained_edges_due_to_rearrangement", "gained_edges_due_to_adjacent_paralogs",
                        "lost_edges_due_to_gene_loss", "lost_edges_due_to_rearrangement"]:
                stratigraphy[key] = 0
                    
            is_retained = dict()
            
            for gene, adjacent_gene, edge_attributes in graph.edges(data=True):
                if gene.parent is None or adjacent_gene is None:
                    stratigraphy["gained_edges_due_to_gene_gain"] += 1
                elif parent_graph.has_edge(gene.parent, adjacent_gene.parent):
                    if parent_occurrence[gene.parent] > 1 and parent_occurrence[adjacent_gene.parent] > 1:
                        stratigraphy["duplicated_edges"] += 1
                    else:
                        stratigraphy["retained_edges"] += 1
                    is_retained[(gene.parent, adjacent_gene.parent)] = True
                    is_retained[(adjacent_gene.parent, gene.parent)] = True
                elif gene.parent == adjacent_gene.parent:
                    stratigraphy["gained_edges_due_to_adjacent_paralogs"] += 1
                else:
                    stratigraphy["gained_edges_due_to_rearrangement"] += 1
            
            for gene, adjacent_gene, edge_attributes in parent_graph.edges(data=True):
                if not (gene, adjacent_gene) in is_retained:
                    try:
                        desc_a = gene.get_at_level(genome)
                        desc_b = adjacent_gene.get_at_level(genome)
                        stratigraphy["lost_edges_due_to_rearrangement"] += 1
                    except:
                        stratigraphy["lost_edges_due_to_gene_loss"] += 1
                    
            tree_node.add_feature('stratigraphy', stratigraphy)      
