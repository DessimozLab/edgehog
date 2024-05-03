#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from . import __version__ as version
import argparse
from .check_args import check_args
from .process_hogs import map_hogs_onto_tree, get_hogxml_entries
from .characterize_species_tree import characterize_tree
from .init_extant_synteny_graphs import init_extant_graphs, init_extant_graphs_from_hdf5
from .infer_ancestral_synteny_graphs import leaves_to_root_synteny_propagation, root_to_leaves_edge_trimming, linearize_graphs
from .write_output import write_output, write_as_hdf5
from .add_ons import date_edges, phylostratify
import time



def main():
    arg_parser = argparse.ArgumentParser(description='edgehog is a software tool that infers an ancestral synteny graph '
                                         'at each internal node of an input species phylogenetic tree')
    arg_parser.add_argument('--version', action='version', help='print version number and exit', version=version)
    arg_parser.add_argument('--output_directory', default='./edgehog_output', type=str, help='path to output directory (default is ./edgehog_output)')
    arg_parser.add_argument('--species_tree', type=str, required=True, help='path to species/genomes phylogenetic tree (newick format)')
    arg_parser.add_argument('--hogs', type=str, required=True, help='path to the HierarchicalGroups.orthoxml file in which HOGs are stored')
    arg_parser.add_argument('--gff_directory', type=str, help='path to directory with the gffs of extant genomes '
                            '(each gff file must be named according to the name of an extant genome / leaf on the species tree)')
    arg_parser.add_argument('--hdf5', type=str, help='path to the hdf5 file (alternative to gff_directory to run edgeHOG on the entire OMA database)')
    arg_parser.add_argument('--date_edges', action='store_true', help='whether the age of edges in extant species should be predicted')
    arg_parser.add_argument('--phylostratify', action='store_true', help='whether the number of edge retention, gain and loss should be analyzed for each node of the species tree')
    arg_parser.add_argument('--max_gaps', type=int, help='max_gaps can be seen as the theoritical maximal number of consecutive novel genes that can emerge between two older genes (default = 3), '
                            'e.g.  if max_gaps = 2: the probabilistic A-b-c-D-E-f-g-h-I-J graph will be turn into A-D-E ; I-J in the ancestor'
                            'while if max_gaps = 3: the probabilistic A-b-c-D-E-f-g-h-I-J graph will be turn into A-D-E-I-J   in the ancestor', default=3)
    # arg_parser.add_argument('--cpu', type=str, default='1', help='number of CPUs to use (default is 1)')
    arg_parser.add_argument('--include_extant_genes', action='store_true', help='include extant genes in output file for ancestral reconstructions.')
    arg_parser.add_argument("--out-format", choices=("TSV", "HDF5"), default="TSV",
                            help="define output format. Can be TSV (tab seperated files) or HDF5 (compatible for integration into oma hdf5)")
    args = arg_parser.parse_args()

    timer = dict()

    start_time = time.time()
    out_dir = check_args(args)

    ham = map_hogs_onto_tree(args.hogs, args.species_tree, args.hdf5)
    hogxml_entries, protein_id_to_hogxml_entry = get_hogxml_entries(ham)
    characterize_tree(ham.taxonomy.tree)
    end_time = time.time()
    timer["preprocessing - loading orthoxml"] = end_time - start_time

    start_time = time.time()
    if args.gff_directory:
        init_extant_graphs(1, ham, args.hogs, args.gff_directory, hogxml_entries, protein_id_to_hogxml_entry)
    elif args.hdf5:
        init_extant_graphs_from_hdf5(ham, args.hdf5)

    end_time = time.time()
    timer["preprocessing - loading extant synteny"] = end_time - start_time

    start_time = time.time()
    leaves_to_root_synteny_propagation(ham, args.max_gaps)
    end_time = time.time()
    timer["bottom-up phase"] = end_time - start_time

    start_time = time.time()
    root_to_leaves_edge_trimming(ham)
    end_time = time.time()
    timer["top-down phase"] = end_time - start_time

    if args.date_edges:
        start_time = time.time()
        date_edges(ham)
        end_time = time.time()
        timer["edge datation"] = end_time - start_time

    start_time = time.time()
    linearize_graphs(ham)
    end_time = time.time()
    timer["linearization"] = end_time - start_time

    if args.phylostratify:
        start_time = time.time()
        phylostratify(ham)
        end_time = time.time()
        timer["phylostratigraphy"] = end_time - start_time

    start_time = time.time()
    if args.out_format == "HDF5":
        write_as_hdf5(args, ham, out_dir)
    else:
        write_output(args, ham, out_dir)
    end_time = time.time()
    timer["writing output"] = end_time - start_time

    print('Completed!')

    print('###################################')
    for step in timer:
        print("%30s:\tdone in %.3f seconds" % (step, timer[step]))


if __name__=='__main__':
    main()

