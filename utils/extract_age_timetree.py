import argparse
import os
from pathlib import Path

import ete3
import pandas as pd


def build_arg_parser():
    """Handle the parameter sent when executing the script from the terminal.

    Returns
    -------
    A parser object with the chosen option and parameters

    """
    parser = argparse.ArgumentParser(description=
            "Make a new version of EdgeHOG extant species result with age of nodes.")
    parser.add_argument("-s", "--species_tree", type=str, help="Original species tree.",
                         required=True)
    parser.add_argument("-t", "--timetree", type=str,
                         help="TimeTree species tree with same namespace (but not necessarily identical leafset)",  # noqa: E501
                         required=True)

    parser.add_argument("-i", "--input_folder",
                        type=str,
                        help="Input folder (EdgeHOG's output).",
                        default="./edgehog_output")
    parser.add_argument("-o", "--output_folder",
                        type=str,
                        help="Output folder (default: corrected_EH).",
                        default="./timed_edgehog/")

    return parser.parse_args()


def read_tree(tree_path : str) -> ete3.Tree:
    return ete3.Tree(tree_path, format=1, quoted_node_names=True)

def uniformize_leafset(tree1 : ete3.Tree , tree2 : ete3.Tree)-> tuple[ete3.Tree,ete3.Tree]:
    l_tree1 = [x.name for x in (tree1.get_leaves())]
    l_tree2 = [x.name for x in (tree2.get_leaves())]
    l_keep = set(l_tree1).intersection(set(l_tree2))
    tree1.prune(l_keep)
    tree2.prune(l_keep)
    return tree1 , tree2

def obtain_age_clade(ori_tree : ete3.Tree, timetree : ete3.Tree, imputation=True) -> dict:
    #Go through the OMA tree and for every clade with clear correspondance, obtain the
    # clade age according to TimeTree (total tre length - branch length of ancestral node to leaves).
    # Give an age of 0 to leaves.
    #In case of disagreement, take the age of the oldest leafset that contain all species in the oma tree
    ori_tree_pruned, timetree = uniformize_leafset(ori_tree.copy(), timetree)
    age_info = {}
    valid_timetree = {}
    leaf= timetree.get_leaves()[0]
    tot_tree_age = leaf.get_distance(timetree)
    for node in ori_tree_pruned.traverse():
        leave_set  = [x.name for x in node.get_leaves()]
        ca = timetree.get_common_ancestor(leave_set)
        if sorted([x.name for x in ca.get_leaves()])==sorted(leave_set):
            age = tot_tree_age - ca.get_distance(timetree)
            age_info[node.name] = float(age)
            valid_timetree[node.name] = True
        elif node.is_leaf():
            age_info[node.name] = 0
        else:
            #Disagreement between taxonomies - for now opt for the longest time
            age = tot_tree_age - ca.get_distance(timetree)
            age_info[node.name] = float(age)
    if imputation:
        age_info = impute_missing(ori_tree, age_info)
    return age_info


def impute_missing(ori_tree : ete3.Tree, age_info : dict):
    #Go back to the full, non-pruned OMA Tree and impute the age of the missing clades
    #as the median between its parent clade age and the age of the oldest children
    for node in ori_tree.traverse("postorder"):
        if node.name not in age_info:
            if node.is_leaf():
                age_info[node.name] = 0

            else:
                ancestor_nr = 1
                cur_anc = node.up
                while cur_anc.name not in age_info:
                    cur_anc = cur_anc.up
                    ancestor_nr+=1
                age_info[node.name] = (age_info[cur_anc.name] + ancestor_nr*max([age_info[x.name] for x in node.children])) / (ancestor_nr+1)
    return age_info

def update_eh_with_ages(extant_file : str, new_extant_file : str, age_info : dict) -> None:
    edgehog_df = pd.read_csv(extant_file, sep="\t")
    edgehog_df["age_MiY"] = [float(age_info[x]) if not pd.isna(x) else pd.NA for x in edgehog_df["predicted_edge_lca"]]
    edgehog_df["age_MiY"] =  edgehog_df["age_MiY"].astype(dtype="Float64")
    edgehog_df.to_csv(new_extant_file, sep="\t")


if __name__ == "__main__":
    args = build_arg_parser()
    #Reading original OMA species tree
    original_tree = args.species_tree
    timetree = args.timetree
    folder = args.input_folder
    newfolder = args.output_folder
    ori_tree = read_tree(original_tree)
    time_tree = read_tree(timetree)
    age_info = obtain_age_clade(ori_tree, time_tree)
    for f in os.listdir(folder):
        if f.endswith("extant_synteny_graph_edges.tsv"):

            path = Path(folder) / f
            new_path = Path(newfolder,f)
            update_eh_with_ages(path, new_path, age_info)
