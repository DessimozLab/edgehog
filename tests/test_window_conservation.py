"""Tests for utils/window_conservation.py and the before_output hook of edgehog.main."""

import importlib.util
import os

import networkx as nx
import pytest

from edgehog import edgehog as edgehog_module
from synteny_helpers import Gene, genes

TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(TESTS)
TEST_DATA = os.path.join(ROOT, "test_data")

_spec = importlib.util.spec_from_file_location("window_conservation",
                                               os.path.join(ROOT, "utils", "window_conservation.py"))
wc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wc)


def path(*nodes):
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    nx.add_path(graph, nodes)
    return graph


def sides(result):
    return sorted([[str(n) for n in side] for side in result])


# --- walk ------------------------------------------------------------------

def test_walk_middle_of_path():
    g = genes("a b c d e")
    graph = path(*g.values())
    assert sides(wc.walk(graph, g["c"], 10)) == [["b", "a"], ["d", "e"]]


def test_walk_limits_to_n():
    g = genes("a b c d e")
    graph = path(*g.values())
    assert sides(wc.walk(graph, g["c"], 1)) == [["b"], ["d"]]


def test_walk_path_end_and_singleton():
    g = genes("a b c")
    graph = path(*g.values())
    assert sides(wc.walk(graph, g["a"], 10)) == [[], ["b", "c"]]
    graph.add_node(Gene("lonely"))
    assert wc.walk(graph, [n for n in graph if n.name == "lonely"][0], 10) == [[], []]


def test_walk_ring_counts_every_node_once():
    g = genes("a b c d e")
    graph = path(*g.values())
    graph.add_edge(g["e"], g["a"])
    left, right = wc.walk(graph, g["a"], 10)
    assert sorted(map(str, left + right)) == ["b", "c", "d", "e"]


def test_walk_skipped_nodes_do_not_use_slots():
    g = genes("a x1 x2 b c")
    graph = path(*g.values())
    result = wc.walk(graph, g["c"], 2, keep=lambda n: not n.name.startswith("x"))
    assert sides(result) == [[], ["b", "a"]]


# --- branch_windows --------------------------------------------------------

def rows(parent_graph, child_graph, window=10):
    return sorted(((str(h), str(c), counts["n_parent_window"], counts["n_child_window"], counts["n_shared"])
                   for h, c, counts in wc.branch_windows(parent_graph, child_graph, window)),
                  key=lambda r: (r[0], r[1]))


def children_of(parents, names):
    """Child genes 'a1'.. whose parent is parents['a'], or None for names starting with 'x'."""
    return {n: Gene(n, None if n.startswith("x") else parents[n[0]]) for n in names.split()}


def test_identical_order_shares_everything():
    p = genes("a b c d")
    c = children_of(p, "a1 b1 c1 d1")
    out = rows(path(*p.values()), path(*c.values()))
    assert all(n_p == n_c == n_s for _, _, n_p, n_c, n_s in out)
    assert [r[:3] for r in out] == [("a", "a1", 3), ("b", "b1", 3), ("c", "c1", 3), ("d", "d1", 3)]


def test_novel_genes_are_skipped_in_child_window():
    p = genes("a b c")
    c = children_of(p, "a1 x1 x2 b1 x3 c1")
    out = rows(path(*p.values()), path(*c.values()), window=1)
    assert ("b", "b1", 2, 2, 2) in out


def test_rearrangement_reduces_shared():
    p = genes("a b c d e")
    c = children_of(p, "a1 b1 c1 e1 d1")   # d-e swapped
    out = {r[0]: r for r in rows(path(*p.values()), path(*c.values()), window=1)}
    assert out["c"] == ("c", "c1", 2, 2, 1)   # parent {b,d}, child {b,e}


def test_duplication_gives_one_row_per_copy_and_excludes_tandem_copy():
    p = genes("a b c")
    c = children_of(p, "a1 b1 b2 c1")       # tandem duplication of b
    out = [r for r in rows(path(*p.values()), path(*c.values()), window=1) if r[0] == "b"]
    # b1's neighbours {a, b2->b}: b itself is excluded, so only {a}
    assert out == [("b", "b1", 2, 1, 1), ("b", "b2", 2, 1, 1)]


def test_loss_gives_row_without_child():
    p = genes("a b c")
    c = children_of(p, "a1 c1")
    out = rows(path(*p.values()), path(*c.values()))
    assert ("b", "None", 2, 0, 0) in out


def test_unmapped_children_are_reported_and_skipped():
    p = genes("a b")
    stray = Gene("stray_hog")
    c = children_of(p, "a1 b1")
    s1 = Gene("s1", stray)
    child = path(c["a1"], s1, c["b1"])
    assert wc.unmapped_children(path(*p.values()), child) == [s1]
    assert ("a", "a1", 1, 1, 1) in rows(path(*p.values()), child)


# --- pipeline hook ---------------------------------------------------------

@pytest.fixture(scope="module")
def test_data_ham(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("edgehog_output")
    seen = {}

    def spy(ham, args, out):
        seen.setdefault("calls", []).append(out)
        seen["files_before"] = os.listdir(out)
        seen["ham"] = ham

    edgehog_module.main(["--hogs", os.path.join(TEST_DATA, "FastOMA_HOGs.orthoxml"),
                         "--species_tree", os.path.join(TEST_DATA, "species_tree.nwk"),
                         "--gff_directory", os.path.join(TEST_DATA, "gff3"),
                         "--output_directory", str(out_dir)], before_output=spy)
    seen["out_dir"] = out_dir
    return seen


def test_hook_called_once_before_output(test_data_ham):
    assert len(test_data_ham["calls"]) == 1
    assert test_data_ham["files_before"] == []
    assert "genome_dict.tsv" in os.listdir(test_data_ham["out_dir"])


def test_window_conservation_on_test_data(test_data_ham):
    df = wc.window_conservation(test_data_ham["ham"], window=10)
    assert list(df.columns) == wc.COLUMNS
    assert len(df) > 0
    assert (df.n_shared <= df[["n_parent_window", "n_child_window"]].min(axis=1)).all()
    assert (df.n_parent_window <= 20).all() and (df.n_child_window <= 20).all()
    # every branch of the tree is covered, including those to leaves
    tree = test_data_ham["ham"].taxonomy.tree
    n_branches = sum(1 for n in tree.traverse() if n.up is not None)
    assert df.groupby(["parent_level", "child_level"]).ngroups == n_branches
    # parent labels match the ones in the linearized output files
    assert df.parent_hog.str.startswith(("HOG_", "rootHOG_")).all()
