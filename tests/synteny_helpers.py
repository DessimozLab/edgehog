"""Lightweight stand-ins for the pyham objects used by the inference phases.

The bottom-up, top-down and linearization phases only need graph nodes with a
``parent`` attribute and tree nodes whose ``genome`` is an instance of
``pyham.AncestralGenome`` / ``pyham.ExtantGenome``. These helpers build such
objects without having to write an orthoXML file.
"""

import types

import networkx as nx
from ete3 import Tree

from edgehog import infer_ancestral_synteny_graphs


class Gene:
    """Graph node with a ``parent`` (the HOG at the parent level, or None)."""

    def __init__(self, name, parent=None):
        self.name = name
        self.parent = parent

    def __repr__(self):
        return self.name


def genes(names, parent=None):
    """Return a dict name -> Gene for whitespace separated names."""
    return {n: Gene(n, parent) for n in names.split()}


class FakeAncestralGenome:
    def __init__(self, name, genes=()):
        self.name = name
        self.genes = list(genes)


class FakeExtantGenome:
    def __init__(self, name):
        self.name = name


def patch_pyham(monkeypatch):
    """Make the isinstance checks in the inference module accept the fake genomes."""
    monkeypatch.setattr(infer_ancestral_synteny_graphs, "pyham",
                        types.SimpleNamespace(AncestralGenome=FakeAncestralGenome,
                                              ExtantGenome=FakeExtantGenome))


def make_ham(newick, ancestral_genes=None):
    """Build a fake Ham around an ete3 tree with named internal nodes.

    Internal nodes get a FakeAncestralGenome (with the genes listed in
    ``ancestral_genes[name]``), leaves get a FakeExtantGenome.
    """
    ancestral_genes = ancestral_genes or {}
    tree = Tree(newick, format=1)
    for node in tree.traverse():
        if node.is_leaf():
            node.add_feature("genome", FakeExtantGenome(node.name))
        else:
            node.add_feature("genome", FakeAncestralGenome(node.name, ancestral_genes.get(node.name, ())))
    tree.n_internal_nodes = sum(1 for n in tree.traverse() if not n.is_leaf())
    nodes = {n.name: n for n in tree.traverse()}
    return types.SimpleNamespace(taxonomy=types.SimpleNamespace(tree=tree)), nodes


def add_edge(graph, a, b, weight=1, u=0, c=0, d=0, **attrs):
    graph.add_edge(a, b, weight=weight, unidirectional=u, convergent=c, divergent=d, **attrs)


def chain_graph(nodes, weights=None, children=None):
    """Linear graph over ``nodes`` with the attributes of an extant synteny graph."""
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    weights = weights or [1] * (len(nodes) - 1)
    for (a, b), w in zip(zip(nodes, nodes[1:]), weights):
        add_edge(graph, a, b, weight=w, children=list(children or [None]), extant_descendants=list(children or [None]))
    return graph


def edge_set(graph):
    """Edges as a set of frozensets of node names (order independent)."""
    return {frozenset((str(a), str(b))) for a, b in graph.edges()}


def edges(*pairs):
    """Build an edge_set-compatible set from strings like 'A-B'."""
    return {frozenset(p.split("-")) for p in pairs}
