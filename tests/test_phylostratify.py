"""Phylostratigraphy (add_ons.phylostratify) on a small synthetic dataset with real pyham.

Species tree (((A,B)AB,(C,D)CD)Root,E)Top. Family p is duplicated on the branch
Top -> Root into the paralogs p1 and p2; p1 is lost on the branch Root -> AB.
Family g is gained in AB (its HOG has no parent at Root).

    A, B   f1 f2 g  p2 f4
    C, D   f1 f2 p1 p2 f4
    E      f1 f2 p  f4

At Root, f2-p1 is kept by the top-down phase because the outgroup E supports it,
so on the branch Root -> AB this edge is lost because p1 is lost.
"""

import pytest
import pyham

from edgehog import add_ons
from edgehog.characterize_species_tree import characterize_tree
from edgehog.infer_ancestral_synteny_graphs import (leaves_to_root_synteny_propagation,
                                                     root_to_leaves_edge_trimming)
from edgehog.init_extant_synteny_graphs import init_extant_graphs
from edgehog.process_hogs import get_hogxml_entries, map_hogs_onto_tree
from gff_helpers import gene, write_gff

SPECIES_TREE = "(((A,B)AB,(C,D)CD)Root,E)Top;"
GENOMES = {
    "A": ["f1", "f2", "g", "p2", "f4"],
    "B": ["f1", "f2", "g", "p2", "f4"],
    "C": ["f1", "f2", "p1", "p2", "f4"],
    "D": ["f1", "f2", "p1", "p2", "f4"],
    "E": ["f1", "f2", "p", "f4"],
}


def write_orthoxml(path):
    ids = {}
    species_xml = []
    for i, (sp, names) in enumerate(sorted(GENOMES.items())):
        genes_xml = []
        for j, name in enumerate(names):
            ids[(sp, name)] = gid = str(1000 * (i + 1) + j)
            genes_xml.append('<gene id="%s" protId="%s_%s"/>' % (gid, sp, name))
        species_xml.append('<species name="%s" NCBITaxId="%d" taxonId="%d"><database name="test" version="1"><genes>%s'
                           '</genes></database></species>' % (sp, i + 1, i + 1, "".join(genes_xml)))

    def group(gid, level, *members):
        return ('<orthologGroup id="%s"><property name="TaxRange" value="%s"/>%s</orthologGroup>'
                % (gid, level, "".join(members)))

    def ref(sp, name):
        return '<geneRef id="%s"/>' % ids[(sp, name)]

    def root_group(gid, name, clades=("AB", "CD")):
        species = {"AB": "AB", "CD": "CD"}
        return group(gid, "Root", *(group("%s_%s" % (gid, c), c, *(ref(sp, name) for sp in species[c]))
                                    for c in clades))

    groups = [group(f, "Top", root_group(f + "_Root", f), ref("E", f)) for f in ("f1", "f2", "f4")]
    groups.append(group("p", "Top", ref("E", "p"),
                        '<paralogGroup>%s%s</paralogGroup>'
                        % (group("p1_Root", "Root", group("p1_CD", "CD", ref("C", "p1"), ref("D", "p1"))),
                           root_group("p2_Root", "p2"))))
    groups.append(group("g", "AB", ref("A", "g"), ref("B", "g")))
    path.write_text('<?xml version="1.0" encoding="utf-8"?>\n'
                    '<orthoXML xmlns="http://orthoXML.org/2011/" origin="test" originVersion="1" version="0.3">'
                    '%s<groups>%s</groups></orthoXML>' % ("".join(species_xml), "".join(groups)))


@pytest.fixture(scope="module")
def ham(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("phylostratify")
    (tmp / "tree.nwk").write_text(SPECIES_TREE)
    write_orthoxml(tmp / "hogs.orthoxml")
    gff_dir = tmp / "gff"
    gff_dir.mkdir()
    for sp, names in GENOMES.items():
        write_gff(gff_dir / ("%s.gff" % sp),
                  *(gene("chr1", 1000 * k + 1, 1000 * k + 500, "+", "%s_%s_gene" % (sp, name), "%s_%s" % (sp, name))
                    for k, name in enumerate(names)))
    ham = map_hogs_onto_tree(str(tmp / "hogs.orthoxml"), str(tmp / "tree.nwk"), None)
    hogxml_entries, protein_id_to_entry = get_hogxml_entries(ham)
    characterize_tree(ham.taxonomy.tree)
    init_extant_graphs(1, ham, None, str(gff_dir), hogxml_entries, protein_id_to_entry, False)
    leaves_to_root_synteny_propagation(ham, 3, False)
    root_to_leaves_edge_trimming(ham, False)
    return ham


@pytest.fixture(scope="module")
def strata(ham):
    add_ons.phylostratify(ham)
    return {n.name: n.props["stratigraphy"] for n in ham.taxonomy.tree.traverse() if n.up}


def family(hog):
    families = {g.prot_id.split("_", 1)[1] for g in hog.get_all_descendant_genes()}
    assert len(families) == 1
    return families.pop()


def family_edges(graph):
    return {frozenset((family(a), family(b))) for a, b in graph.edges()}


def test_dataset_has_edge_lost_with_paralog(ham):
    nodes = {n.name: n for n in ham.taxonomy.tree.traverse()}
    root_edges = family_edges(nodes["Root"].props["top_down_synteny"])
    assert frozenset(("f2", "p1")) in root_edges
    assert {frozenset(("f2", "g")), frozenset(("g", "p2"))} <= family_edges(nodes["AB"].props["top_down_synteny"])


def test_edge_lost_with_paralog_counts_as_gene_loss(ham, strata):
    # f2-p1 (and p1-p2, if kept at Root) disappear in AB because p1 is lost there.
    # The paralog p2 survives in AB, which must not make the loss look like a rearrangement.
    nodes = {n.name: n for n in ham.taxonomy.tree.traverse()}
    root_edges = family_edges(nodes["Root"].props["top_down_synteny"])
    with_p1 = sum(1 for e in root_edges if "p1" in e)
    assert strata["AB"]["lost_edges_due_to_gene_loss"] == with_p1
    # f2-p2 is broken up by the insertion of g: both genes are still present in AB
    assert strata["AB"]["lost_edges_due_to_rearrangement"] == 1


def test_edges_to_gained_gene_count_as_gene_gain(strata):
    # f2-g and g-p2, whichever endpoint networkx yields first
    assert strata["AB"]["gained_edges_due_to_gene_gain"] == 2
    assert strata["AB"]["gained_edges_due_to_rearrangement"] == 0


def test_retained_edges(strata):
    assert strata["AB"]["retained_edges"] == 2  # f1-f2, p2-f4
    assert strata["A"]["retained_edges"] == 4
    for key in ("lost_edges_due_to_gene_loss", "lost_edges_due_to_rearrangement", "gained_edges_due_to_gene_gain"):
        assert strata["A"][key] == 0


def test_does_not_search_whole_families(ham, monkeypatch):
    # HOG.get_at_level walks the whole top-level family per call, which made
    # phylostratify take days on large datasets.
    calls = []
    monkeypatch.setattr(pyham.abstractgene.AbstractGene, "get_at_level", lambda *args: calls.append(args))
    add_ons.phylostratify(ham)
    assert calls == []
