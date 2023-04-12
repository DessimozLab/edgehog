#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 11:59:21 2022

@author: charles
"""

def compute_sum_dist_to_desc_leaves(node, sum_d = 0, n_leaves = 0, n_internal_nodes = 0):
    n_internal_nodes = 0
    if node.is_leaf():
        n_leaves = 1
        sum_d = node.get_distance(node.up)
        node.n_desc_leaves = 0
        node.sum_dist_to_desc_leaves = 0
        node.n_internal_nodes = 0
    else:
        n_leaves = 0
        sum_d = 0
        for child in node.children:
            res = compute_sum_dist_to_desc_leaves(child, sum_d, n_leaves, n_internal_nodes)
            sum_d += res[0] 
            n_leaves += res[1]
            n_internal_nodes += res[2]
        n_internal_nodes += 1
        node.sum_dist_to_desc_leaves = sum_d
        node.n_desc_leaves = n_leaves
        node.n_internal_nodes = n_internal_nodes
        if node.up:
            sum_d += node.get_distance(node.up) * n_leaves
    if node.up:
        return sum_d, n_leaves, n_internal_nodes
      

def compute_red_score(node, red = 0, level_from_root = 0):
    if node.is_leaf():
        red = 1
    elif not node.up:
        red = 0
    else:
        p = red
        d = node.get_distance(node.up)
        u = (node.sum_dist_to_desc_leaves + (d * node.n_desc_leaves)) / node.n_desc_leaves
        red = p + (d/u) * (1-p)
    node.level = level_from_root
    node.red = red
    for child in node.children:
        compute_red_score(child, red, level_from_root + 1)
        

# the list of nodes per level will be useful for subsequent parallelization of the code
# def get_nodes_per_level_from_root(node, level_from_root, nodes_per_level):
#     node.level = level_from_root
#     if not node.is_leaf():
#         if not level_from_root in nodes_per_level:
#             nodes_per_level[level_from_root] = list()
#         nodes_per_level[level_from_root].append(node)
#         children = node.get_children()
#         for child in children:
#             get_nodes_per_level_from_root(child, level_from_root + 1, nodes_per_level)


def characterize_tree(tree):
    compute_sum_dist_to_desc_leaves(tree)
    compute_red_score(tree)
    
    
