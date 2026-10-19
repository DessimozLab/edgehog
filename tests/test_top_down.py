"""Tests for the top-down phase (root_to_leaves_edge_trimming).

Reference: Supplementary Information of Bernard et al. 2025, Nat Ecol Evol,
section "Parsimonious ancestral synteny graphs": an edge is kept if it is
present in the propagated network of at least two children, or present in the
propagated network of one child and its propagated form is present in the
pruned network of the parent.
"""

import networkx as nx
import pytest

from edgehog.infer_ancestral_synteny_graphs import root_to_leaves_edge_trimming
from synteny_helpers import Gene, add_edge, edge_set, edges, make_ham, patch_pyham


def _graph(nodes, *edge_specs):
    """edge_specs: (a, b, weight[, u, c, d])"""
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    for spec in edge_specs:
        a, b, w = spec[:3]
        u, c, d = spec[3:] if len(spec) > 3 else (0, 0, 0)
        add_edge(graph, a, b, weight=w, u=u, c=c, d=d)
    return graph


# ---------------------------------------------------------------------------
# single ancestral level: R with leaves L1, L2 (and L3 for polytomies)
# ---------------------------------------------------------------------------

@pytest.fixture
def root_only(monkeypatch):
    patch_pyham(monkeypatch)

    def build(newick, bottom_up, transitive):
        ham, nodes = make_ham(newick)
        nodes["R"].add_prop("bottom_up_synteny", bottom_up)
        for leaf, graph in transitive.items():
            nodes[leaf].add_prop("transitive_synteny", graph)
        return ham, nodes
    return build


@pytest.mark.parametrize("orient", [False, True])
def test_edge_supported_by_two_children_is_kept(root_only, orient):
    P = {i: Gene("P%d" % i) for i in range(1, 4)}
    ham, nodes = root_only("(L1,L2)R;",
                           _graph(P.values(), (P[1], P[2], 2, 2, 0, 0), (P[2], P[3], 1, 1, 0, 0)),
                           {"L1": _graph(P.values(), (P[1], P[2], 1, 1, 0, 0), (P[2], P[3], 1, 1, 0, 0)),
                            "L2": _graph(P.values(), (P[1], P[2], 1, 1, 0, 0))})
    root_to_leaves_edge_trimming(ham, orient)
    graph = nodes["R"].props["top_down_synteny"]
    # P2-P3 only in L1 and there is no outgroup at the root -> removed
    assert edge_set(graph) == edges("P1-P2")
    assert set(graph.nodes) == set(P.values())
    # bottom-up graph is left untouched
    assert edge_set(nodes["R"].props["bottom_up_synteny"]) == edges("P1-P2", "P2-P3")


def test_polytomy_edge_supported_by_two_of_three_children(root_only):
    P = {i: Gene("P%d" % i) for i in range(1, 3)}
    ham, nodes = root_only("(L1,L2,L3)R;",
                           _graph(P.values(), (P[1], P[2], 2)),
                           {"L1": _graph(P.values(), (P[1], P[2], 1)),
                            "L2": _graph(P.values()),
                            "L3": _graph(P.values(), (P[1], P[2], 1))})
    root_to_leaves_edge_trimming(ham, False)
    assert edge_set(nodes["R"].props["top_down_synteny"]) == edges("P1-P2")


def test_edge_with_weight_from_single_child_is_removed(root_only):
    # weight 2 but both counts come from the same child (e.g. tandem duplicates)
    P = {i: Gene("P%d" % i) for i in range(1, 3)}
    ham, nodes = root_only("(L1,L2)R;",
                           _graph(P.values(), (P[1], P[2], 2)),
                           {"L1": _graph(P.values(), (P[1], P[2], 2)),
                            "L2": _graph(P.values())})
    root_to_leaves_edge_trimming(ham, False)
    assert nodes["R"].props["top_down_synteny"].number_of_edges() == 0


# Issue #1. Paper: "When the --orient_edges option is enabled, the trimming
# algorithm applies to each edge context. Only parsimoniously supported contexts
# are retained, with unsupported ones reset to 0. If an edge is supported but
# none of its contexts are [...], context weights remain unchanged."
# Code: u, c, d are computed but never written back to the graph
# (infer_ancestral_synteny_graphs.py:372-378).

@pytest.mark.xfail(strict=True, reason="#1: unsupported orientation contexts are not reset in the top-down phase")
def test_orientation_unsupported_context_is_reset(root_only):
    P = {i: Gene("P%d" % i) for i in range(1, 3)}
    ham, nodes = root_only("(L1,L2)R;",
                           _graph(P.values(), (P[1], P[2], 3, 2, 1, 0)),
                           {"L1": _graph(P.values(), (P[1], P[2], 1, 1, 0, 0)),
                            "L2": _graph(P.values(), (P[1], P[2], 2, 1, 1, 0))})
    root_to_leaves_edge_trimming(ham, True)
    data = nodes["R"].props["top_down_synteny"][P[1]][P[2]]
    # co-directional seen in both children, convergent only in L2
    assert (data["unidirectional"], data["convergent"], data["divergent"]) == (2, 0, 0)


def test_orientation_unchanged_if_no_context_supported(root_only):
    P = {i: Gene("P%d" % i) for i in range(1, 3)}
    ham, nodes = root_only("(L1,L2)R;",
                           _graph(P.values(), (P[1], P[2], 2, 0, 1, 1)),
                           {"L1": _graph(P.values(), (P[1], P[2], 1, 0, 1, 0)),
                            "L2": _graph(P.values(), (P[1], P[2], 1, 0, 0, 1))})
    root_to_leaves_edge_trimming(ham, True)
    data = nodes["R"].props["top_down_synteny"][P[1]][P[2]]
    assert (data["unidirectional"], data["convergent"], data["divergent"]) == (0, 1, 1)


# ---------------------------------------------------------------------------
# two ancestral levels: ((L1,L2)N,L3)R -- support from the parent (outgroup)
# ---------------------------------------------------------------------------

@pytest.fixture
def two_levels(monkeypatch):
    """Ancestral genes Q1..Q7 at R and P1..P7 (+P7b) at N, Pi.parent == Qi.

    P5 has no parent (gained in N); P7 and P7b are in-paralogs from Q7.
    In N, only L1 supports the edges except P6-P7 (L1 and L2).
    """
    patch_pyham(monkeypatch)
    ham, nodes = make_ham("((L1,L2)N,L3)R;")
    Q = {i: Gene("Q%d" % i) for i in range(1, 8)}
    P = {i: Gene("P%d" % i, Q[i]) for i in range(1, 8)}
    P[5].parent = None
    P["7b"] = Gene("P7b", Q[7])

    nodes["R"].add_prop("bottom_up_synteny", _graph(
        Q.values(), (Q[1], Q[2], 2), (Q[2], Q[3], 1), (Q[4], Q[6], 2)))
    nodes["N"].add_prop("transitive_synteny", _graph(
        Q.values(), (Q[1], Q[2], 1), (Q[2], Q[3], 1), (Q[4], Q[6], 1)))
    nodes["L3"].add_prop("transitive_synteny", _graph(
        Q.values(), (Q[1], Q[2], 1), (Q[4], Q[6], 1)))

    single = [(P[1], P[2], 1), (P[2], P[3], 1), (P[4], P[5], 1), (P[5], P[6], 1), (P[7], P["7b"], 1)]
    nodes["N"].add_prop("bottom_up_synteny", _graph(P.values(), *single, (P[6], P[7], 2)))
    nodes["N"].add_prop("real_edge_to_transitive_edge", {
        (P[4], P[5]): (Q[4], Q[6]), (P[5], P[4]): (Q[4], Q[6]),
        (P[5], P[6]): (Q[4], Q[6]), (P[6], P[5]): (Q[4], Q[6])})
    nodes["L1"].add_prop("transitive_synteny", _graph(P.values(), *single, (P[6], P[7], 1)))
    nodes["L2"].add_prop("transitive_synteny", _graph(P.values(), (P[6], P[7], 1)))
    return ham, nodes, P, Q


def test_root_level_of_two_levels(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, False)
    assert edge_set(nodes["R"].props["top_down_synteny"]) == edges("Q1-Q2", "Q4-Q6")


def test_single_child_edge_supported_by_parent_is_kept(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, False)
    graph = nodes["N"].props["top_down_synteny"]
    assert graph.has_edge(P[1], P[2])        # Q1-Q2 in pruned parent graph
    assert not graph.has_edge(P[2], P[3])    # Q2-Q3 was pruned at R
    assert graph.has_edge(P[6], P[7])        # two children


def test_single_child_edge_supported_through_transitive_edge(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, False)
    graph = nodes["N"].props["top_down_synteny"]
    # P4-P5-P6 was propagated as Q4-Q6 (P5 gained in N), which is kept at R
    assert graph.has_edge(P[4], P[5]) and graph.has_edge(P[5], P[6])


def test_adjacent_in_paralogs_supported_by_parental_gene(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, False)
    assert nodes["N"].props["top_down_synteny"].has_edge(P[7], P["7b"])


def test_nb_internal_nodes_without_orientation(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, False)
    graph = nodes["N"].props["top_down_synteny"]
    assert graph[P[6]][P[7]]["nb_internal_nodes_from_ancestor_with_updated_weight"] == 0
    assert graph[P[1]][P[2]]["nb_internal_nodes_from_ancestor_with_updated_weight"] == 1


def test_same_edges_kept_with_and_without_orientation(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, False)
    without = {n: edge_set(nodes[n].props["top_down_synteny"]) for n in ("R", "N")}
    root_to_leaves_edge_trimming(ham, True)
    assert without == {n: edge_set(nodes[n].props["top_down_synteny"]) for n in ("R", "N")}


# Issue #5: with --orient_edges the children loop does not break and the parent
# check always runs, which overwrites the value 0 set for edges supported by two
# children (infer_ancestral_synteny_graphs.py:338-346).

@pytest.mark.xfail(strict=True, reason="#5: nb_internal_nodes_from_ancestor_with_updated_weight depends on --orient_edges")
def test_nb_internal_nodes_independent_of_orientation(two_levels):
    ham, nodes, P, Q = two_levels
    root_to_leaves_edge_trimming(ham, True)
    graph = nodes["N"].props["top_down_synteny"]
    assert graph[P[6]][P[7]]["nb_internal_nodes_from_ancestor_with_updated_weight"] == 0
