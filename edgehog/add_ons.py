#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pyham


def date_edges(ham):
    print('###################################')
    print('Predicting age of extant edges ...')
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        genome = tree_node.props.get('genome')
        if genome is None:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if isinstance(genome, pyham.AncestralGenome):
            graph = tree_node.props['top_down_synteny']
        elif isinstance(genome, pyham.ExtantGenome):
            graph = tree_node.props['bottom_up_synteny']
        if tree_node.up:
            parent_graph = tree_node.up.props['top_down_synteny']
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
            graph[gene][adjacent_gene]["lca"] = graph[gene][adjacent_gene]["lca_node"].props['genome'].name
            graph[gene][adjacent_gene]["age"] = 1 - graph[gene][adjacent_gene]["lca_node"].props['red']


def phylostratify(ham):
    print('###################################')
    print('Phylostratigraphy ...')
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        genome = tree_node.props.get('genome')
        if genome is None:
            print('No genome stored for {}'.format(tree_node.name))
            continue
        if tree_node.up:
            if isinstance(genome, pyham.AncestralGenome):
                graph = tree_node.props['top_down_synteny']
            elif isinstance(genome, pyham.ExtantGenome):
                graph = tree_node.props['bottom_up_synteny']
            parent_graph = tree_node.up.props['top_down_synteny']
            
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
                if gene.parent is None or adjacent_gene.parent is None:
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
            
            # A parent HOG survives on this branch if it has a direct child at this
            # level (pyham inserts the intermediate HOGs, so no deeper search is
            # needed). HOG.get_at_level would instead search the whole top-level
            # family, which is slow and also finds the descendants of paralogs.
            has_descendant = dict()

            def survives(hog):
                if hog not in has_descendant:
                    has_descendant[hog] = any(child.genome == genome for child in hog.children)
                return has_descendant[hog]

            for gene, adjacent_gene in parent_graph.edges():
                if not (gene, adjacent_gene) in is_retained:
                    if survives(gene) and survives(adjacent_gene):
                        stratigraphy["lost_edges_due_to_rearrangement"] += 1
                    else:
                        stratigraphy["lost_edges_due_to_gene_loss"] += 1
                    
            tree_node.add_prop('stratigraphy', stratigraphy)      
