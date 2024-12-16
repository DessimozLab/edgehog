import ete3
import os
import argparse
import pandas as pd

def build_arg_parser():
    """Handle the parameter sent when executing the script from the terminal

    Returns
    -----------
    A parser object with the chosen option and parameters"""

    parser = argparse.ArgumentParser(description="Make a new version of EdgeHOG extant species result with age of nodes.")
    parser.add_argument('-s', '--species_tree', type=str, help="Original species tree.", required=True)    
    parser.add_argument('-t', '--timetree', type=str,  help='TimeTree species tree with same namespace (but not necessarily identical leafset)', required=True)
    
    parser.add_argument('-i', '--input_folder', type=str, help="Input folder (EdgeHOG's output).", default='./edgehog_output')
    parser.add_argument('-o', '--output_folder', type=str, help="Output folder (default: corrected_EH).", default="./timed_edgehog/")

    args = parser.parse_args()
    return args


def read_tree(tree_path):
    tree =ete3.Tree(tree_path, format=1, quoted_node_names=True)
    return tree

def uniformize_leafset(tree1, tree2):
    l_tree1 = [x.name for x in (tree1.get_leaves())]
    l_tree2 = [x.name for x in (tree2.get_leaves())]
    l_keep = set(l_tree1).intersection(set(l_tree2))
    tree1.prune(l_keep)
    tree2.prune(l_keep)
    return tree1 , tree2

def obtain_age_clade(ori_tree, timetree, imputation=True):
    #Go through the OMA tree and for every clade with clear correspondance, obtain the clade age according to TimeTree (total tre length - branch length of ancestral node to leaves). Give a age of 0 to leaves.
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


def impute_missing(ori_tree, age_info):
    #Go back to the full, non-pruned OMA Tree and impute the age of the missing clades as the median between its parent clade age and the age of the oldest children
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

def update_EH_with_ages(extant_file, new_extant_file, age_info):
    df = pd.read_csv(extant_file, sep='\t')
    df['age_MiY'] = [float(age_info[x]) if not pd.isna(x) else pd.NA for x in df['predicted_edge_lca']]
    df['age_MiY'] =  df['age_MiY'].astype(dtype='Float64')
    df.to_csv(new_extant_file, sep="\t")


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
            path = os.path.join(folder,f)
            new_path = os.path.join(newfolder,f)
            update_EH_with_ages(path, new_path, age_info)