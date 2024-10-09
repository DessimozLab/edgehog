#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy
import numpy as np
import pandas as pd
import os
import networkx as nx
import pyham
import tables


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


class HDF5Writer:
    def __init__(self, fname, oma_db_fn, date_edges=False, orient_edges=False):
        self.fname = fname
        self.oma_db = oma_db_fn
        self.date_edges = date_edges
        self.orient_edges = orient_edges

    def __enter__(self):
        self.oma_h5 = tables.open_file(self.oma_db, 'r')
        self.h5 = tables.open_file(self.fname, 'w', filters=tables.Filters(complevel=6, complib="blosc"))
        self.tax2taxid = self._load_taxname_2_taxid()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.h5.close()
        self.oma_h5.close()

    def _load_taxname_2_taxid(self):
        tab = self.oma_h5.get_node("/Taxonomy")
        return {row['Name'].decode(): int(row['NCBITaxonId']) for row in tab.read()}

    def _load_hog_at_level(self, taxid):
        try:
            tab = self.oma_h5.get_node('/AncestralGenomes/tax{}/Hogs'.format(taxid))
        except tables.NoSuchNodeError:
            tab = self.oma_h5.get_node('/Hogs_per_Level/tax{}'.format(taxid))
        return {row['ID'].decode(): row_nr for row_nr, row in enumerate(tab.read())}

    def add_graph_at_level(self, taxid, tree_node):
        try:
            from pyoma.browser.tablefmt import AncestralSyntenyRels
        except ImportError:
            print(f"pyoma library is required to write output as HDF5 files. "
                  f"Please install edgehog with the `oma` extra activated, i.e. `pip install edgehog[oma]`.")
            import sys
            sys.exit(2)
        hogid_lookup = self._load_hog_at_level(taxid)
        dtype = tables.dtype_from_descr(AncestralSyntenyRels)
        orient_enum = AncestralSyntenyRels.columns['Orientation'].enum
        dfs = []
        for evidence, graph in zip(
                ("linearized", "parsimonious", "any"),
                (tree_node.linear_synteny, tree_node.top_down_synteny, tree_node.bottom_up_synteny)):
            data = []
            ev = AncestralSyntenyRels.columns['Evidence'].enum[evidence]
            print(f"process level {taxid} - graph {evidence} - |N|,|V| = {len(graph.nodes)},{len(graph.edges)}")

            for u, v, edge_data in graph.edges.data():
                w = edge_data.get("weight", default=1)
                h1 = u.hog_id.rsplit('_')[0]
                h2 = v.hog_id.rsplit('_')[0]
                lca, orient, orient_score = -1, -1, 0
                if self.date_edges:
                    try:
                        lca_clade = edge_data['lca']
                        lca = self.tax2taxid[lca_clade]
                    except KeyError:
                        pass
                if self.orient_edges:
                    keys = ['unidirectional', 'divergent', 'convergent']
                    try:
                        orient_weights = numpy.array(list(map(lambda o: edge_data[o], keys)))
                    except KeyError:
                        pass
                    else:
                        orient = orient_enum[keys[numpy.argmax(orient_weights)]]
                        orient_score = orient_weights.max()
                rec = (hogid_lookup[h1], hogid_lookup[h2], w, ev, lca, orient, orient_score)
                data.append(rec)
            df = pd.DataFrame.from_records(numpy.array(data, dtype=dtype))
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        sumdf = df.groupby(by=["HogRow1", "HogRow2"], as_index=False).agg({
            "Weight": "max",
            "Evidence": "min",
            # -1 indicates 'n/a', all other values should be consistent.
            "LCA_taxid": lambda x: x[x != -1].max() if not x[x != -1].empty else -1,
            "Orientation": lambda x: x[x != -1].max() if not x[x != -1].empty else -1,
            "OrientationScore": "max",

        })
        as_array = sumdf.to_records(index=False, column_dtypes={"LCA_taxid": np.int32,})
        tab = self.h5.create_table("/AncestralGenomes/tax{}".format(taxid),
                                   "Synteny", description=AncestralSyntenyRels,
                                   obj=as_array,
                                   expectedrows=len(as_array),
                                   createparents=True)
        for col in ("HogRow1", "HogRow2", "Weight", "Evidence", "LCA_taxid", ):
            tab.colinstances[col].create_csindex()
        print(f" hdf5 level {taxid} written: {len(tab)} rows.")

    def get_taxid_from_hog_names(self, graph):
        taxids = set([])
        for u in graph.nodes:
            try:
                hogid, taxid = u.hog_id.rsplit('_')
            except ValueError:
                continue
            taxids.add(taxid)
        if len(taxids) > 1:
            raise Exception("multiple taxids on this level: {}".format(taxids))
        elif len(taxids) == 0:
            return 0
        return taxids.pop()


def write_as_hdf5(args, ham, out_dir):
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

