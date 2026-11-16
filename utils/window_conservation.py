#!/usr/bin/env python3
"""Conservation of gene neighbourhoods along the branches of the species tree.

For every branch parent -> child and every HOG ``h`` in the parent's linearized
synteny graph, the window of up to ``W`` genes on each side of ``h`` is compared
with the window around each descendant ``c`` of ``h`` in the child's graph
(``linear_synteny`` for ancestral nodes, the extant graph for leaves). Child
genes without a parent HOG (gained on the branch) are skipped when building the
child window, and the remaining child genes are represented by their parent
HOG, so that both windows are sets of parent-level HOGs.

One row is written per (h, c) pair. If ``h`` has several descendants at the child
level (duplication on the branch), it gets one row per copy. If it has none (loss),
it gets one row with an empty ``child_hog``. Only counts are reported:

    n_parent_window   distinct HOGs in the parent window of h
    n_child_window    distinct parent HOGs in the window of c (h itself excluded)
    n_shared          size of the intersection of both windows
    n_*_left/right    the same per side (the side orientation is arbitrary, so
                      left/right of parent and child are not necessarily aligned)

Usage: run edgehog as usual, with the extra options of this script in front:

    uv run python utils/window_conservation.py --window 10 \\
        --analysis-out window_conservation.tsv.gz <edgehog arguments>

The edgehog output files are still written as usual.
"""

import argparse
import sys

import pandas as pd
import pyham

from edgehog.edgehog import main as edgehog_main
from edgehog.write_output import label_nodes

COLUMNS = ["parent_level", "child_level", "parent_hog", "child_hog",
           "n_parent_window", "n_child_window", "n_shared",
           "n_parent_left", "n_parent_right", "n_child_left", "n_child_right"]


def walk(graph, start, n, keep=None):
    """Return the first ``n`` nodes kept on each side of ``start`` in a graph of max degree 2.

    Nodes for which ``keep(node)`` is false are passed over without counting.
    A walk stops at the end of a path, or when it comes back to ``start`` or to a
    node already collected (ring components). The result has two lists, one per
    neighbour of ``start``. The second list is empty if ``start`` has fewer than two neighbours.
    """
    sides = []
    seen = {start}
    for first in graph.neighbors(start):
        side = []
        prev, node = start, first
        while node is not None and node not in seen and len(side) < n:
            seen.add(node)
            if keep is None or keep(node):
                side.append(node)
            nxt = None
            for nb in graph.neighbors(node):
                if nb != prev:
                    nxt = nb
                    break
            prev, node = node, nxt
        sides.append(side)
    if len(sides) > 2:
        raise ValueError("node {!r} has degree {} in a linear graph".format(start, len(sides)))
    while len(sides) < 2:
        sides.append([])
    return sides


def branch_windows(parent_graph, child_graph, window, parent_of=lambda g: g.parent):
    """Yield ``(h, c, counts)`` for one branch, ``c`` is None if ``h`` has no descendant.

    Child nodes whose parent HOG is not in ``parent_graph`` are treated like
    novel genes (skipped); ``unmapped_children`` lists them.
    """
    children_of = {}
    for c in child_graph.nodes:
        p = parent_of(c)
        if p is not None and p in parent_graph:
            children_of.setdefault(p, []).append(c)

    def mapped(x):
        p = parent_of(x)
        return p is not None and p in parent_graph

    for h in parent_graph.nodes:
        p_left, p_right = (set(side) for side in walk(parent_graph, h, window))
        p_all = p_left | p_right
        children = children_of.get(h)
        if not children:
            yield h, None, dict(n_parent_window=len(p_all), n_child_window=0, n_shared=0,
                                n_parent_left=len(p_left), n_parent_right=len(p_right),
                                n_child_left=0, n_child_right=0)
            continue
        for c in children:
            c_left, c_right = ({parent_of(x) for x in side} - {h}
                               for side in walk(child_graph, c, window, keep=mapped))
            c_all = c_left | c_right
            yield h, c, dict(n_parent_window=len(p_all), n_child_window=len(c_all),
                             n_shared=len(p_all & c_all),
                             n_parent_left=len(p_left), n_parent_right=len(p_right),
                             n_child_left=len(c_left), n_child_right=len(c_right))


def unmapped_children(parent_graph, child_graph, parent_of=lambda g: g.parent):
    """Child nodes with a parent HOG that is not a node of ``parent_graph``."""
    return [c for c in child_graph.nodes
            if parent_of(c) is not None and parent_of(c) not in parent_graph]


def linear_graph(tree_node):
    genome = tree_node.props.get('genome')
    if genome is None:
        return None
    if isinstance(genome, pyham.AncestralGenome):
        return tree_node.props['linear_synteny']
    return tree_node.props['bottom_up_synteny']


def labels(tree_node):
    """Gene labels as used in the edgehog TSV outputs of this node."""
    genome = tree_node.props['genome']
    if isinstance(genome, pyham.AncestralGenome):
        # write_output labels the ancestral graphs from the bottom-up graph
        label_dict, _ = label_nodes(tree_node.props['bottom_up_synteny'])
        return label_dict.__getitem__
    return lambda gene: gene.prot_id


def window_conservation(ham, window=10):
    """DataFrame with one row per (parent HOG, child HOG) over all branches of the tree."""
    rows = []
    label_cache = {}

    def label_of(tree_node):
        if tree_node not in label_cache:
            label_cache[tree_node] = labels(tree_node)
        return label_cache[tree_node]

    for child_node in ham.taxonomy.tree.traverse('preorder'):
        parent_node = child_node.up
        if parent_node is None:
            continue
        parent_graph, child_graph = linear_graph(parent_node), linear_graph(child_node)
        if parent_graph is None or child_graph is None:
            continue
        parent_name = parent_node.props['genome'].name
        child_name = child_node.props['genome'].name
        unmapped = unmapped_children(parent_graph, child_graph)
        if unmapped:
            print("window_conservation: {} genes of {} have a parent HOG that is not in the graph of {}, "
                  "they are ignored".format(len(unmapped), child_name, parent_name), file=sys.stderr)
        parent_label, child_label = label_of(parent_node), label_of(child_node)
        for h, c, counts in branch_windows(parent_graph, child_graph, window):
            row = dict(parent_level=parent_name, child_level=child_name,
                       parent_hog=parent_label(h), child_hog=None if c is None else child_label(c))
            row.update(counts)
            rows.append(row)
    return pd.DataFrame(rows, columns=COLUMNS)


def run(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], add_help=False)
    parser.add_argument('--window', type=int, default=10,
                        help='number of genes on each side of a HOG (default 10)')
    parser.add_argument('--analysis-out', required=True,
                        help='path of the TSV file to write (gzipped if it ends in .gz)')
    opts, edgehog_argv = parser.parse_known_args(argv)

    def analyse(ham, args, out_dir):
        print('###################################')
        print('Window conservation analysis (window = {}) ...'.format(opts.window))
        df = window_conservation(ham, opts.window)
        df.to_csv(opts.analysis_out, sep='\t', index=False)
        print('wrote {} rows to {}'.format(len(df), opts.analysis_out))

    edgehog_main(edgehog_argv, before_output=analyse)


if __name__ == '__main__':
    run()
