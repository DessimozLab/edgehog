"""Tests for the linearization phase (linearize_graphs).

Reference: Supplementary Information of Bernard et al. 2025, Nat Ecol Evol,
section "Linearized graphs", and Fig. 2c of the main text.
"""

import random

import networkx as nx
import pytest

from edgehog.infer_ancestral_synteny_graphs import compute_neighbor_cost, linearize_graphs
from synteny_helpers import edge_set, edges, make_ham, patch_pyham


@pytest.fixture
def linearize(monkeypatch):
    patch_pyham(monkeypatch)

    def run(graph):
        ham, nodes = make_ham("(L1,L2)R;")
        nodes["R"].add_feature("top_down_synteny", graph)
        linearize_graphs(ham)
        assert nodes["R"].top_down_synteny is graph  # input must not be modified
        return nodes["R"].linear_synteny
    return run


def weighted(*specs):
    """Graph from strings like 'A-B:3'."""
    graph = nx.Graph()
    for spec in specs:
        pair, w = spec.split(":")
        a, b = pair.split("-")
        graph.add_edge(a, b, weight=int(w))
    return graph


def test_compute_neighbor_cost_follows_linear_path():
    graph = weighted("X-B:2", "B-B2:1", "B2-B3:4", "B3-B4:1", "B3-B5:1")
    # path X->B->B2->B3 stops at B3 (degree 3)
    assert compute_neighbor_cost(graph, "B", "X") == 2 + 1 + 4
    assert compute_neighbor_cost(graph, "B4", "B3") == 1


def test_compute_neighbor_cost_terminates_on_cycle():
    graph = weighted("A-B:1", "B-C:1", "C-A:1")
    assert compute_neighbor_cost(graph, "B", "A", max_path_len=10) == 10


def test_linear_graph_is_unchanged(linearize):
    graph = weighted("A-B:1", "B-C:5", "D-E:2")
    graph.add_node("F")
    result = linearize(graph)
    assert edge_set(result) == edge_set(graph)
    assert set(result.nodes) == set(graph.nodes)


def test_paper_example_neighbor_priority(linearize):
    # Supplementary Information: conflicting node with neighbors A (weight 3),
    # B (weight 2, along a path to another node) and terminal C (weight 2);
    # the priority list for edge removal is C, B, A.
    graph = weighted("X-A:3", "X-B:2", "B-B2:1", "X-C:2")
    assert edge_set(linearize(graph)) == edges("X-A", "X-B", "B-B2")


def test_lowest_weight_edges_removed_first(linearize):
    graph = weighted("H-A:5", "H-B:4", "H-C:1", "H-D:2", "H-E:3")
    assert edge_set(linearize(graph)) == edges("H-A", "H-B")


def test_reintroduction_of_non_reciprocal_removals(linearize):
    # example from the code comment in linearize_graphs:
    #     A---4---B---4---C                         A---4---B---4---C
    #         2 /   \ 2          -Linearization->
    # W---4---X---2---Y---4---Z                   W---4---X---2---Y---4---Z
    graph = weighted("A-B:4", "B-C:4", "B-X:2", "B-Y:2", "W-X:4", "X-Y:2", "Y-Z:4")
    assert edge_set(linearize(graph)) == edges("A-B", "B-C", "W-X", "X-Y", "Y-Z")


# Issue #3. Paper: removed edges are reintroduced by priority of their weight
# (criterion 1), then of the weight of the resulting linear connected component
# (criterion 2). Code: np.argsort(weights)[::-1] reverses the preceding
# cost ordering for equal weights (infer_ancestral_synteny_graphs.py:498-502).

@pytest.mark.xfail(strict=True, reason="#3: tie-break on component weight is inverted when reintroducing edges")
def test_reintroduction_prefers_heavier_component(linearize):
    # Y and Z each drop their weight-2 edge to X; hub B then drops X, Y and Z.
    # Both X-Y and X-Z can be reintroduced (weight 2), but only one of them.
    # X-Z yields the heavier component (W-X-Z-Z2: 3+2+9) than X-Y (W-X-Y-Y2: 3+2+5).
    graph = weighted("B-A:9", "B-C:9", "B-X:3", "B-Y:3", "B-Z:3",
                     "X-W:3", "X-Y:2", "X-Z:2", "Y-Y2:5", "Z-Z2:9")
    result = linearize(graph)
    assert result.has_edge("X", "Z")
    assert not result.has_edge("X", "Y")


@pytest.mark.parametrize("seed", range(20))
def test_random_graphs_become_linear_subgraphs(linearize, seed):
    rng = random.Random(seed)
    graph = nx.gnm_random_graph(40, 70, seed=seed)
    for a, b in graph.edges():
        graph[a][b]["weight"] = rng.randint(1, 5)
    result = linearize(graph)
    assert max(d for _, d in result.degree()) <= 2
    assert set(result.edges()) <= set(graph.edges())
    assert set(result.nodes) == set(graph.nodes)
