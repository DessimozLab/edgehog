"""Tests for the bottom-up phase (infer_ancestral_synteny_graphs.py).

Reference: Supplementary Information of Bernard et al. 2025, Nat Ecol Evol,
section "Ancestral synteny graphs".
"""

import itertools

import networkx as nx
import pytest

from edgehog.infer_ancestral_synteny_graphs import (build_transitive_graph,
                                                     leaves_to_root_synteny_propagation,
                                                     orient_reconnected_edges)
from synteny_helpers import Gene, add_edge, edge_set, edges, make_ham, patch_pyham


# ---------------------------------------------------------------------------
# orient_reconnected_edges: merging the contexts of A-b and b-C into A-C
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ctx_a, ctx_b, expected", [
    # i) two co-directional contexts -> co-directional
    ("u", "u", "u"),
    # ii) divergent + convergent (or vice versa) -> co-directional
    ("d", "c", "u"),
    ("c", "d", "u"),
    # iii) divergent + co-directional (or vice versa) -> divergent
    ("d", "u", "d"),
    ("u", "d", "d"),
    # iv) convergent + co-directional (or vice versa) -> convergent
    ("c", "u", "c"),
    ("u", "c", "c"),
])
def test_orient_reconnected_edges_rules(ctx_a, ctx_b, expected):
    def ctx(name, w):
        return [w if k == name else 0 for k in "ucd"]
    u, c, d = orient_reconnected_edges(3, *ctx(ctx_a, 3), 2, *ctx(ctx_b, 2))
    result = dict(zip("ucd", (u, c, d)))
    assert result[expected] == pytest.approx(2)  # min(wA, wB)
    assert sum(result.values()) == pytest.approx(2)


@pytest.mark.parametrize("ctx", ["c", "d"])
def test_orient_reconnected_edges_incompatible_contexts(ctx):
    # -> <- followed by -> <- (or <- -> then <- ->) is impossible for the middle gene
    args = [1 if k == ctx else 0 for k in "ucd"]
    assert orient_reconnected_edges(1, *args, 1, *args) == (0, 0, 0)


def test_orient_reconnected_edges_mixed_contexts_scaled_to_min_weight():
    u, c, d = orient_reconnected_edges(3, 2, 1, 0, 2, 2, 0, 0)
    assert u > 0 and c > 0 and d == 0
    assert u + c + d == pytest.approx(2)


# ---------------------------------------------------------------------------
# build_transitive_graph: removing genes without parent in a child graph
# ---------------------------------------------------------------------------

def _transitive_input(names, weights=None, extra_edges=()):
    """Chain over names; upper-case genes get a parent HOG, lower-case ones not."""
    if isinstance(names, str):
        names = names.split()
    gs ={n: Gene(n, Gene("P" + n) if n.isupper() else None) for n in names}
    graph = nx.Graph()
    graph.add_nodes_from(gs[n] for n in names)
    weights = weights or [1] * (len(names) - 1)
    for i, (a, b) in enumerate(zip(names, names[1:])):
        add_edge(graph, gs[a], gs[b], weight=weights[i], u=weights[i],
                 children=["child%d" % i], extant_descendants=["leaf%d" % i])
    for a, b, w in extra_edges:
        add_edge(graph, gs[a], gs[b], weight=w, u=w, children=["extra"], extant_descendants=["extra"])
    return graph, gs


def _with_node_order(graph, order):
    """Same graph, but with nodes inserted in the given order (controls iteration order)."""
    g = nx.Graph()
    g.add_nodes_from(order)
    g.add_edges_from(graph.edges(data=True))
    return g


def test_transitive_reconnects_gene_without_parent():
    graph, g = _transitive_input("A b C", weights=[3, 2])
    build_transitive_graph(graph, orient_edges=False, max_gap_length=2)
    assert g["b"] not in graph
    assert edge_set(graph) == edges("A-C")
    data = graph[g["A"]][g["C"]]
    assert data["weight"] == 2  # minimal weight of the two removed edges
    assert sorted(data["children"]) == ["child0", "child1"]
    assert sorted(data["extant_descendants"]) == ["leaf0", "leaf1"]


def test_transitive_keeps_edges_between_genes_with_parent():
    graph, g = _transitive_input("A B C", weights=[3, 2])
    build_transitive_graph(graph, orient_edges=False, max_gap_length=2)
    assert edge_set(graph) == edges("A-B", "B-C")
    assert graph[g["A"]][g["B"]]["weight"] == 3


def test_transitive_terminal_gene_without_parent_is_not_bridged():
    # the gene stays in the working copy, but its edge has no ancestral counterpart
    graph, g = _transitive_input("A b")
    mapping = build_transitive_graph(graph, orient_edges=False, max_gap_length=2)
    assert edge_set(graph) == edges("A-b")
    assert mapping == {}


def test_transitive_no_reconnection_around_branching_gene():
    # only parentless genes of degree 2 are bridged; a degree-3 gene is a break point
    gs = {n: Gene(n, None if n == "b" else Gene("P" + n)) for n in "AbCD"}
    graph = nx.Graph()
    for other in "ACD":
        add_edge(graph, gs["b"], gs[other], children=[None], extant_descendants=[None])
    mapping = build_transitive_graph(graph, orient_edges=False, max_gap_length=2)
    assert edge_set(graph) == edges("b-A", "b-C", "b-D")
    assert mapping == {}


@pytest.mark.parametrize("names, max_gap, expected", [
    ("A b c D", 2, edges("A-D")),
    ("A b c d E", 2, set()),
    ("A b c d E", 3, edges("A-E")),
    ("A b c D E f g h I J", 2, edges("A-D", "D-E", "I-J")),
    ("A b c D E f g h I J", 3, edges("A-D", "D-E", "E-I", "I-J")),
])
def test_transitive_max_gap_independent_of_node_order(names, max_gap, expected):
    names = names.split()
    graph, g = _transitive_input(names)
    parentless = [n for n in names if n.islower()]
    for perm in itertools.permutations(parentless):
        order = [g[n] for n in perm] + [g[n] for n in names if n.isupper()]
        work = _with_node_order(graph, order)
        build_transitive_graph(work, orient_edges=False, max_gap_length=max_gap)
        assert edge_set(work) == expected, "node order %s" % (perm,)


@pytest.mark.parametrize("order", [("b", "c"), ("c", "b")])
def test_transitive_maps_all_real_edges_of_a_gap(order):
    graph, g = _transitive_input("A b c D".split())
    work = _with_node_order(graph, [g[n] for n in order] + [g["A"], g["D"]])
    mapping = build_transitive_graph(work, orient_edges=False, max_gap_length=2)
    expected = {g["A"].parent, g["D"].parent}
    for a, b in [("A", "b"), ("b", "c"), ("c", "D")]:
        assert set(mapping[(g[a], g[b])]) == expected
        assert set(mapping[(g[b], g[a])]) == expected


def test_transitive_reconnected_edge_orientation():
    graph, g = _transitive_input("A b C", weights=[3, 2])  # both edges co-directional
    build_transitive_graph(graph, orient_edges=True, max_gap_length=2)
    data = graph[g["A"]][g["C"]]
    assert (data["unidirectional"], data["convergent"], data["divergent"]) == (2, 0, 0)


# Issue #2. Paper: "in case the reconnected edge was already present in the
# network, this minimal weight is added to the initial edge weight".
# Code: build_transitive_graph only reconnects if not graph.has_edge(...),
# see infer_ancestral_synteny_graphs.py:135.

@pytest.mark.xfail(strict=True, reason="#2: weight of an already existing reconnected edge is not incremented")
def test_transitive_existing_edge_weight_is_incremented():
    graph, g = _transitive_input("A b C", weights=[3, 2], extra_edges=[("A", "C", 5)])
    build_transitive_graph(graph, orient_edges=False, max_gap_length=2)
    assert graph[g["A"]][g["C"]]["weight"] == 5 + 2


@pytest.mark.xfail(strict=True, reason="#2: real edges are not mapped when the reconnected edge already exists")
def test_transitive_existing_edge_real_edges_are_mapped():
    graph, g = _transitive_input("A b C", weights=[3, 2], extra_edges=[("A", "C", 5)])
    mapping = build_transitive_graph(graph, orient_edges=False, max_gap_length=2)
    assert set(mapping[(g["A"], g["b"])]) == {g["A"].parent, g["C"].parent}


# ---------------------------------------------------------------------------
# leaves_to_root_synteny_propagation
# ---------------------------------------------------------------------------

@pytest.fixture
def two_leaf_tree(monkeypatch):
    """Root R with leaves L1 and L2; ancestral genes P1..P5.

    L1: a1 - a2 - x - a3     (x has no parent HOG: gene gained in L1)
        a2 - a2b - a1        (a2b is an in-paralog of a2, i.e. also maps to P2)
    L2: b1 - b2 ; b3 - b4    (two contigs)
    P5 has no descendant edges and must remain a singleton.
    """
    patch_pyham(monkeypatch)
    P = {i: Gene("P%d" % i) for i in range(1, 6)}
    ham, nodes = make_ham("(L1,L2)R;", {"R": list(P.values())})
    a = {n: Gene(n, P[int(n[1])] if n != "x" else None) for n in ["a1", "a2", "a2b", "x", "a3"]}
    b = {n: Gene(n, P[int(n[1])]) for n in ["b1", "b2", "b3", "b4"]}

    g1 = nx.Graph()
    add_edge(g1, a["a1"], a["a2"], u=1, children=[None], extant_descendants=[None])
    add_edge(g1, a["a2"], a["x"], u=1, children=[None], extant_descendants=[None])
    add_edge(g1, a["x"], a["a3"], c=1, children=[None], extant_descendants=[None])
    add_edge(g1, a["a2b"], a["a1"], d=1, children=[None], extant_descendants=[None])
    add_edge(g1, a["a2"], a["a2b"], u=1, children=[None], extant_descendants=[None])
    g2 = nx.Graph()
    add_edge(g2, b["b1"], b["b2"], u=1, children=[None], extant_descendants=[None])
    add_edge(g2, b["b3"], b["b4"], u=1, children=[None], extant_descendants=[None])
    nodes["L1"].add_prop("bottom_up_synteny", g1)
    nodes["L2"].add_prop("bottom_up_synteny", g2)
    return ham, nodes, P


def test_propagation_sums_weights_over_children(two_leaf_tree):
    ham, nodes, P = two_leaf_tree
    leaves_to_root_synteny_propagation(ham, max_gaps=2, orient_edges=True)
    graph = nodes["R"].props["bottom_up_synteny"]
    assert edge_set(graph) == edges("P1-P2", "P2-P3", "P3-P4")
    p12 = graph[P[1]][P[2]]
    # a1-a2 and a2b-a1 in L1 both map onto P1-P2, plus b1-b2 in L2
    assert p12["weight"] == 3
    assert sorted(c.name for c in p12["children"]) == ["L1", "L1", "L2"]
    assert sorted(c.name for c in p12["extant_descendants"]) == ["L1", "L1", "L2"]
    assert (p12["unidirectional"], p12["convergent"], p12["divergent"]) == (2, 0, 1)


def test_propagation_bridges_gained_gene(two_leaf_tree):
    ham, nodes, P = two_leaf_tree
    leaves_to_root_synteny_propagation(ham, max_gaps=2, orient_edges=True)
    p23 = nodes["R"].props["bottom_up_synteny"][P[2]][P[3]]
    assert p23["weight"] == 1
    # a2 -x-> a3 is co-directional followed by convergent -> convergent
    assert (p23["unidirectional"], p23["convergent"], p23["divergent"]) == (0, 1, 0)


def test_propagation_keeps_singleton_ancestral_genes(two_leaf_tree):
    ham, nodes, P = two_leaf_tree
    leaves_to_root_synteny_propagation(ham, max_gaps=2, orient_edges=False)
    graph = nodes["R"].props["bottom_up_synteny"]
    assert P[5] in graph and graph.degree(P[5]) == 0
    assert not any(a == b for a, b in graph.edges())


def test_propagation_stores_child_transitive_graphs(two_leaf_tree):
    ham, nodes, P = two_leaf_tree
    leaves_to_root_synteny_propagation(ham, max_gaps=2, orient_edges=False)
    t1 = nodes["L1"].props["transitive_synteny"]
    t2 = nodes["L2"].props["transitive_synteny"]
    assert edge_set(t1) == edges("P1-P2", "P2-P3")
    assert t1[P[1]][P[2]]["weight"] == 2
    assert edge_set(t2) == edges("P1-P2", "P3-P4")
    # the per-child weights add up to the parent weight
    parent = nodes["R"].props["bottom_up_synteny"]
    for a, b, data in parent.edges(data=True):
        contributed = sum(t[a][b]["weight"] for t in (t1, t2) if t.has_edge(a, b))
        assert contributed == data["weight"]
