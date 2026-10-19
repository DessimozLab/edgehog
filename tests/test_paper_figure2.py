"""The worked examples of Fig. 2 in Bernard et al. 2025, Nat Ecol Evol
(https://doi.org/10.1038/s41559-025-02818-0).

Fig. 2a/b: species tree ((S1,S2)Pink,S3)Root with HOGs A-E. A duplicated
into A1/A2 between Root and Pink, B emerged in Pink, D and E duplicated in
the S3 lineage, X and Y are gained genes (not in any HOG).

    S1: A2 A1 B C D E
    S2: X A2 E A1 B
    S3: A Y C D1 E1 D2 E2

Fig. 2c: linearization of a synteny network with conflicting nodes Y, W, I, V.
"""

import networkx as nx
import pytest

from edgehog.characterize_species_tree import characterize_tree
from edgehog.infer_ancestral_synteny_graphs import (leaves_to_root_synteny_propagation,
                                                     linearize_graphs,
                                                     root_to_leaves_edge_trimming)
from edgehog.init_extant_synteny_graphs import init_extant_graphs
from edgehog.process_hogs import get_hogxml_entries, map_hogs_onto_tree
from gff_helpers import gene, write_gff
from synteny_helpers import edges, make_ham, patch_pyham

SPECIES_TREE = "((S1,S2)Pink,S3)Root;"
GENOMES = {
    "S1": ["A2", "A1", "B", "C", "D", "E"],
    "S2": ["X", "A2", "E", "A1", "B"],
    "S3": ["A", "Y", "C", "D1", "E1", "D2", "E2"],
}

# HOGs as nested orthoXML groups; {S1_A1} is replaced by a geneRef
GROUPS = """
<orthologGroup id="A"><property name="TaxRange" value="Root"/>
  <paralogGroup>
    <orthologGroup id="A1"><property name="TaxRange" value="Pink"/>{S1_A1}{S2_A1}</orthologGroup>
    <orthologGroup id="A2"><property name="TaxRange" value="Pink"/>{S1_A2}{S2_A2}</orthologGroup>
  </paralogGroup>
  {S3_A}
</orthologGroup>
<orthologGroup id="B"><property name="TaxRange" value="Pink"/>{S1_B}{S2_B}</orthologGroup>
<orthologGroup id="C"><property name="TaxRange" value="Root"/>
  <orthologGroup id="C_Pink"><property name="TaxRange" value="Pink"/>{S1_C}</orthologGroup>
  {S3_C}
</orthologGroup>
<orthologGroup id="D"><property name="TaxRange" value="Root"/>
  <orthologGroup id="D_Pink"><property name="TaxRange" value="Pink"/>{S1_D}</orthologGroup>
  <paralogGroup>{S3_D1}{S3_D2}</paralogGroup>
</orthologGroup>
<orthologGroup id="E"><property name="TaxRange" value="Root"/>
  <orthologGroup id="E_Pink"><property name="TaxRange" value="Pink"/>{S1_E}{S2_E}</orthologGroup>
  <paralogGroup>{S3_E1}{S3_E2}</paralogGroup>
</orthologGroup>
"""


def write_orthoxml(path):
    refs = {}
    species_xml = []
    for i, (sp, names) in enumerate(sorted(GENOMES.items())):
        genes_xml = []
        for j, name in enumerate(names):
            gid = str(1000 * (i + 1) + j)
            refs["%s_%s" % (sp, name)] = '<geneRef id="%s"/>' % gid
            genes_xml.append('<gene id="%s" protId="%s_%s"/>' % (gid, sp, name))
        species_xml.append('<species name="%s" NCBITaxId="%d"><database name="test" version="1"><genes>%s'
                           '</genes></database></species>' % (sp, i + 1, "".join(genes_xml)))
    path.write_text('<?xml version="1.0" encoding="utf-8"?>\n'
                    '<orthoXML xmlns="http://orthoXML.org/2011/" origin="test" originVersion="1" version="0.3">'
                    '%s<groups>%s</groups></orthoXML>' % ("".join(species_xml), GROUPS.format(**refs)))


@pytest.fixture(scope="module")
def fig2(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("fig2")
    (tmp / "tree.nwk").write_text(SPECIES_TREE)
    write_orthoxml(tmp / "hogs.orthoxml")
    gff_dir = tmp / "gff"
    gff_dir.mkdir()
    for sp, names in GENOMES.items():
        write_gff(gff_dir / ("%s.gff" % sp),
                  *[gene("chr1", 1000 * k + 1, 1000 * k + 500, "+", "%s_%s_gene" % (sp, name), "%s_%s" % (sp, name))
                    for k, name in enumerate(names)])

    ham = map_hogs_onto_tree(str(tmp / "hogs.orthoxml"), str(tmp / "tree.nwk"), None)
    hogxml_entries, protein_id_to_entry = get_hogxml_entries(ham)
    characterize_tree(ham.taxonomy.tree)
    init_extant_graphs(1, ham, None, str(gff_dir), hogxml_entries, protein_id_to_entry, False)
    leaves_to_root_synteny_propagation(ham, 3, False)
    root_to_leaves_edge_trimming(ham, False)
    linearize_graphs(ham)
    return {node.name: node for node in ham.taxonomy.tree.traverse()}


def label(hog):
    """Name of an ancestral gene as in the figure (A1/A2 at Pink, A at Root)."""
    names = {g.prot_id.split("_", 1)[1] for g in hog.get_all_descendant_genes()}
    if len(names) > 1:
        names = {n.rstrip("0123456789") for n in names}
    assert len(names) == 1, names
    return names.pop()


def weighted_edges(graph):
    return {frozenset((label(a), label(b))): data["weight"] for a, b, data in graph.edges(data=True)}


# ---------------------------------------------------------------------------
# Fig. 2a: bottom-up phase
# ---------------------------------------------------------------------------

def test_fig2a_pink_bottom_up_graph(fig2):
    # S1 and S2 propagated graphs: A2-A1-B-C-D-E and A2-E-A1-B (X is gained)
    assert weighted_edges(fig2["Pink"].bottom_up_synteny) == {
        frozenset("ED"): 1, frozenset("DC"): 1, frozenset("CB"): 1, frozenset(("B", "A1")): 2,
        frozenset(("A1", "A2")): 1, frozenset(("A2", "E")): 1, frozenset(("E", "A1")): 1}


def test_fig2a_pink_propagated_graph(fig2):
    # propagated forms: A1-B -> A-C, B-C -> A-C, A2-A1 -> A
    assert weighted_edges(fig2["Pink"].transitive_synteny) == {
        frozenset("ED"): 1, frozenset("DC"): 1, frozenset("AC"): 1, frozenset("AE"): 2}


def test_fig2a_pink_propagated_forms(fig2):
    mapping = fig2["Pink"].real_edge_to_transitive_edge
    forms = {frozenset((label(a), label(b))): frozenset(label(n) for n in form) for (a, b), form in mapping.items()}
    assert forms[frozenset(("A1", "B"))] == frozenset("AC")
    assert forms[frozenset(("B", "C"))] == frozenset("AC")


def test_fig2a_s3_propagated_graph_edges(fig2):
    # S3 propagated graph: A-C-D-E (Y is gained, D1/E1/D2/E2 collapse onto D-E)
    assert set(weighted_edges(fig2["S3"].transitive_synteny)) == edges("A-C", "C-D", "D-E")


def test_fig2a_root_bottom_up_graph_edges(fig2):
    assert set(weighted_edges(fig2["Root"].bottom_up_synteny)) == edges("E-D", "D-C", "A-C", "A-E")
    weights = weighted_edges(fig2["Root"].bottom_up_synteny)
    assert weights[frozenset("DC")] == 2
    assert weights[frozenset("AC")] == 2
    assert weights[frozenset("AE")] == 2


# Fig. 2a shows weight 2 for D-E in the S3 propagated graph (and hence 3 at the
# root). The code counts each of the three extant adjacencies D1-E1, E1-D2 and
# D2-E2 (weight 3, and 4 at the root).

@pytest.mark.xfail(strict=True, reason="Fig. 2a: D-E weight from the S3 tandem duplicates is 3 in the code, 2 in the figure")
def test_fig2a_tandem_duplicate_weights(fig2):
    assert weighted_edges(fig2["S3"].transitive_synteny)[frozenset("DE")] == 2
    assert weighted_edges(fig2["Root"].bottom_up_synteny)[frozenset("DE")] == 3


# ---------------------------------------------------------------------------
# Fig. 2b: top-down phase
# ---------------------------------------------------------------------------

def test_fig2b_root_top_down(fig2):
    # E-D, D-C, A-C supported by two children; A-E only by Pink -> removed
    assert set(weighted_edges(fig2["Root"].top_down_synteny)) == edges("E-D", "D-C", "A-C")


def test_fig2b_pink_top_down(fig2):
    # kept: B-A1 (two children), E-D and D-C (one child and parent),
    # C-B (propagated form A-C in parent), A2-A1 (propagated form A in parent);
    # removed: A2-E and E-A1 (propagated form A-E not in parent)
    assert set(weighted_edges(fig2["Pink"].top_down_synteny)) == edges("E-D", "D-C", "C-B", "B-A1", "A1-A2")


def test_fig2b_linearized(fig2):
    assert set(weighted_edges(fig2["Root"].linear_synteny)) == edges("E-D", "D-C", "A-C")
    assert set(weighted_edges(fig2["Pink"].linear_synteny)) == edges("A2-A1", "A1-B", "B-C", "C-D", "D-E")


# ---------------------------------------------------------------------------
# Fig. 2c: linearization
# ---------------------------------------------------------------------------

FIG2C = ["Z-Y:3", "Y-Q:1", "Y-X:2", "Q-W:1", "X-W:2", "W-V:2",
         "U-V:3", "V-G:1", "V-I:2", "I-K:1", "I-J:2", "K-J:3"]


@pytest.fixture
def fig2c(monkeypatch):
    patch_pyham(monkeypatch)
    graph = nx.Graph()
    for spec in FIG2C:
        pair, w = spec.split(":")
        graph.add_edge(*pair.split("-"), weight=int(w))
    ham, nodes = make_ham("(L1,L2)R;")
    nodes["R"].add_feature("top_down_synteny", graph)
    linearize_graphs(ham)
    return {frozenset(e) for e in nodes["R"].linear_synteny.edges()}


def test_fig2c_conflicts_resolved(fig2c):
    # Y drops Q (e=1); W drops Q (e=1); I drops K (e=1);
    # V drops G (e=1) and then I (e=2, path weight 7 < 9 for W)
    assert edges("Z-Y", "Y-X", "X-W", "W-V", "V-U", "I-J", "J-K") <= fig2c
    assert not fig2c & edges("Y-Q", "Q-W", "V-G", "V-I")


# Issue #10 (found with this example): after V is resolved, I and K both have
# a single neighbor, so the earlier removed edge I-K is reintroduced
# (infer_ancestral_synteny_graphs.py:487-512) and I-J-K becomes a cycle.
# The figure shows the linear contig I-J-K.

@pytest.mark.xfail(strict=True, reason="#10: reintroducing a removed edge can close a cycle (I-J-K)")
def test_fig2c_final_contigs(fig2c):
    assert fig2c == edges("Z-Y", "Y-X", "X-W", "W-V", "V-U", "I-J", "J-K")
