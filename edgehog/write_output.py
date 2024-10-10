#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy
import pandas as pd
import os
import networkx as nx
import pyham


def label_nodes(graph):
    label_dict = dict()
    annotation_dict = dict()
    visited_node = dict()
    visited_label = dict()
    subhogs = dict()
    for node in graph.nodes:
        if not node in visited_node:
            visited_node[node] = True
            if node.hog_id:
                label = 'HOG_' + node.hog_id
            else:
                label = 'HOG_' + node.get_top_level_hog().hog_id
            if not node.parent:
                label = 'root' + label
            if label not in visited_label:
                visited_label[label] = True
                subhogs[label] = 0
            else:
                subhogs[label] += 1
                label += "_c" + str(subhogs[label])
            label_dict[node] = label
            #annotation_dict[node] = ";".join([gene.genome.name + ":" + gene.prot_id for gene in node.get_all_descendant_genes()])
            annotation_dict[node] = ";".join([gene.genome.name + ":" + gene.prot_id.split(' ')[0] for gene in node.get_all_descendant_genes()])
    return label_dict, annotation_dict


def graph_to_df(graph, genome, edge_datation, orient_edges, label_dict=None, annotation_dict=None, include_extant_genes=False):
    edge_dict = dict()
    singleton_dict = dict()
    column_names = ['gene1', 'gene2', 'weight', 'contiguous_region',
                    'nb_internal_nodes_from_ancestor_with_updated_weight', 'supporting_children']
    if orient_edges:
        column_names += ['predicted_transcriptional_orientation', 'orientation_score']
    if include_extant_genes:
        column_names += ['gene1_extant_annotations', 'gene2_extant_annotations']
    if edge_datation:
        column_names += ['predicted_edge_age_relative_to_root', 'predicted_edge_lca']
    for c in column_names:
        edge_dict[c] = list()
        singleton_dict[c] = list()
    car_counter = 0
    for cc in nx.connected_components(graph):
        if len(cc) > 1:
            for u, v, edge_data in graph.subgraph(cc).edges.data():
                w = edge_data['weight']
                edge_dict['contiguous_region'].append(car_counter)
                edge_dict['weight'].append(w)
                if isinstance(genome, pyham.AncestralGenome):
                    edge_dict['gene1'].append(label_dict[u])
                    edge_dict['gene2'].append(label_dict[v])

                    if orient_edges:
                        # in case of equalities, here is the hierarchy: unidirectional > divergent > convergent
                        e_u, e_d, e_c = edge_data['unidirectional'], edge_data['divergent'], edge_data['convergent']
                        edge_dict['predicted_transcriptional_orientation'].append(['unidirectional','divergent','convergent',][numpy.argmax([e_u, e_d, e_c])])
                        edge_dict['orientation_score'].append(round(max(e_u,e_c,e_d)))

                    if include_extant_genes:
                        edge_dict['gene1_extant_annotations'].append(annotation_dict[u])
                        edge_dict['gene2_extant_annotations'].append(annotation_dict[v])

                    if 'nb_internal_nodes_from_ancestor_with_updated_weight' in graph[u][v]:
                        edge_dict['nb_internal_nodes_from_ancestor_with_updated_weight'].append(graph[u][v]['nb_internal_nodes_from_ancestor_with_updated_weight'])
                    else:
                        edge_dict['nb_internal_nodes_from_ancestor_with_updated_weight'].append(None)
                    edge_dict['supporting_children'].append(";".join(sorted(set([child.genome.name for child in graph[u][v]['children']]))))
                else:
                    edge_dict['nb_internal_nodes_from_ancestor_with_updated_weight'].append(0)
                    edge_dict['supporting_children'].append(None)
                    edge_dict['gene1'].append(u.prot_id)
                    edge_dict['gene2'].append(v.prot_id)
                    if orient_edges:
                        e_u, e_d, e_c = edge_data['unidirectional'], edge_data['divergent'], edge_data['convergent']
                        edge_dict['predicted_transcriptional_orientation'].append(['unidirectional','divergent','convergent',][numpy.argmax([e_u, e_d, e_c])])
                        edge_dict['orientation_score'].append(round(max(e_u,e_d,e_c)))
                    if include_extant_genes:
                        edge_dict['gene1_extant_annotations'].append(u.prot_id)
                        edge_dict['gene2_extant_annotations'].append(v.prot_id)

                if edge_datation:
                    edge_dict['predicted_edge_lca'].append(graph[u][v]["lca"])
                    edge_dict['predicted_edge_age_relative_to_root'].append('%.2f' % graph[u][v]["age"])
            car_counter += 1
        else:
            u = list(cc)[0]
            for key in ['contiguous_region', 'weight', 'gene2', 'nb_internal_nodes_from_ancestor_with_updated_weight', 'supporting_children']:
                singleton_dict[key].append(None)
            if orient_edges:
                singleton_dict['predicted_transcriptional_orientation'].append(None)
                singleton_dict['orientation_score'].append(None)
            if edge_datation:
                singleton_dict['predicted_edge_lca'].append(None)
                singleton_dict['predicted_edge_age_relative_to_root'].append(None)
            if isinstance(genome, pyham.AncestralGenome):
                singleton_dict['gene1'].append(label_dict[u])
                if include_extant_genes:
                    singleton_dict['gene1_extant_annotations'].append(annotation_dict[u])
                    singleton_dict['gene2_extant_annotations'].append(None)
            else:
                singleton_dict['gene1'].append(u.prot_id)
                if include_extant_genes:
                    singleton_dict['gene1_extant_annotations'].append(u.prot_id)
                    singleton_dict['gene2_extant_annotations'].append(None)

    for c in column_names:
        edge_dict[c] += singleton_dict[c]
    df = pd.DataFrame.from_records(edge_dict)
    return df[column_names]


def write_as_hdf5(args, ham, out_dir):
    from .hdf5 import HDF5Writer
    with HDF5Writer(os.path.join(out_dir, "Synteny.h5"), args.hdf5, date_edges=args.date_edges, orient_edges=args.orient_edges) as writer:
        for tree_node in ham.taxonomy.tree.traverse("preorder"):
            try:
                genome = tree_node.genome
            except AttributeError:
                print('No genome stored for {}'.format(tree_node.name))
                continue
            if isinstance(genome, pyham.AncestralGenome):
                taxid = writer.get_taxid_from_hog_names(tree_node.linear_synteny)
                writer.add_graph_at_level(taxid, tree_node)


def write_output(args, ham, out_dir):
    print('###################################')
    print('Writing output files ...')
    genome_dict = {'genome_id': [], 'name': [], 'nb_descendant_leaves': [], 'level_from_root': [], 'RED_score': []}
    if args.phylostratify:
        stratigraphy_columns = ["retained_edges", "duplicated_edges", "gained_edges_due_to_gene_gain",
                                "gained_edges_due_to_rearrangement", "gained_edges_due_to_adjacent_paralogs",
                                "lost_edges_due_to_gene_loss", "lost_edges_due_to_rearrangement"]
        for key in stratigraphy_columns:
            genome_dict[key] = [None]
    g = 0
    for tree_node in ham.taxonomy.tree.traverse('preorder'):
        try:
            genome = tree_node.genome
        except AttributeError:
            print('No genome stored for {}'.format(tree_node.name))
            continue

        if isinstance(genome, pyham.AncestralGenome):
            label_dict, annotation_dict = label_nodes(tree_node.bottom_up_synteny)

            df = graph_to_df(tree_node.bottom_up_synteny, genome, False, args.orient_edges, label_dict,
                             annotation_dict, args.include_extant_genes)
            df.to_csv(os.path.join(out_dir, str(g) + '_bottom-up_synteny_graph_edges.tsv.gz'),
                      compression='gzip',
                      index=False,
                      header=True,
                      sep='\t')

            df = graph_to_df(tree_node.top_down_synteny, genome, args.date_edges, args.orient_edges, label_dict,
                             annotation_dict, args.include_extant_genes)
            df.to_csv(os.path.join(out_dir, str(g) + '_top-down_synteny_graph_edges.tsv.gz'),
                      compression='gzip',
                      index=False,
                      header=True,
                      sep='\t')

            df = graph_to_df(tree_node.linear_synteny, genome, args.date_edges, args.orient_edges, label_dict,
                             annotation_dict, args.include_extant_genes)
            df.to_csv(os.path.join(out_dir, str(g) + '_linearized_synteny_graph_edges.tsv.gz'),
                      compression='gzip',
                      index=False,
                      header=True,
                      sep='\t')
        else:
            # extant genome
            df = graph_to_df(tree_node.bottom_up_synteny, genome, args.date_edges, args.orient_edges)
            df.to_csv(os.path.join(out_dir, str(g) + '_extant_synteny_graph_edges.tsv.gz'),
                      compression='gzip',
                      index=False,
                      header=True,
                      sep='\t')


        genome_dict['genome_id'].append(g)
        genome_dict['name'].append(genome.name)
        genome_dict['nb_descendant_leaves'].append(tree_node.n_desc_leaves)
        genome_dict['level_from_root'].append(tree_node.level)
        genome_dict['RED_score'].append('%.2f' % tree_node.red)

        if args.phylostratify and tree_node.up:
            stratigraphy = tree_node.stratigraphy
            for key in stratigraphy_columns:
                genome_dict[key].append(stratigraphy[key])

        g += 1

    genome_df = pd.DataFrame.from_records(genome_dict)
    genome_df[['genome_id', 'nb_descendant_leaves', 'level_from_root', 'RED_score', 'name']].to_csv(os.path.join(out_dir, 'genome_dict.tsv'), index = False, header = True, sep = '\t')

    if args.phylostratify:
        genome_df[['genome_id', 'nb_descendant_leaves', 'level_from_root', 'RED_score', 'name'] + stratigraphy_columns].to_csv(os.path.join(out_dir, 'phylostratigraphy.tsv'), index = False, header = True, sep = '\t')

