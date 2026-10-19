"""End-to-end test on a tiny synthetic dataset with known ancestral gene orders.

Species tree ((A,B)AB,(C,D)CD)Root; families f1..f7 are present in all species,
gene x is gained in A (not part of any HOG).

    A    chr1: f1 x f2 f3 f4 f5 f6 f7
    B    chr1: f1 f2 f3 f4 f5 f6 f7
    C    chr1: f1 f2 f3 f4 f5 f6        chr2: f7
    D    chr1: f1 f2 f4 f3 f5 f6        chr2: f7     (inversion of f3-f4)

Expected: f1-...-f6 at every ancestor, f6-f7 only in AB, and the novel
adjacencies of D (f2-f4, f3-f5) nowhere.
"""

import pytest

from edgehog.characterize_species_tree import characterize_tree
from edgehog.infer_ancestral_synteny_graphs import (leaves_to_root_synteny_propagation,
                                                     linearize_graphs,
                                                     root_to_leaves_edge_trimming)
from edgehog.init_extant_synteny_graphs import init_extant_graphs
from edgehog.process_hogs import get_hogxml_entries, map_hogs_onto_tree
from gff_helpers import gene, write_gff
from synteny_helpers import edges

SPECIES_TREE = "((A,B)AB,(C,D)CD)Root;"
GENOMES = {
    "A": [["f1", "x", "f2", "f3", "f4", "f5", "f6", "f7"]],
    "B": [["f1", "f2", "f3", "f4", "f5", "f6", "f7"]],
    "C": [["f1", "f2", "f3", "f4", "f5", "f6"], ["f7"]],
    "D": [["f1", "f2", "f4", "f3", "f5", "f6"], ["f7"]],
}
FAMILIES = ["f%d" % i for i in range(1, 8)]
# SPECIES_TREE as orthoXML <taxonomy>; the leaf ids are the NCBITaxIds of the species
TAXONOMY = ('<taxonomy><taxon id="10" name="Root">'
            '<taxon id="11" name="AB"><taxon id="1" name="A"/><taxon id="2" name="B"/></taxon>'
            '<taxon id="12" name="CD"><taxon id="3" name="C"/><taxon id="4" name="D"/></taxon>'
            '</taxon></taxonomy>')
CHAIN = ("f1-f2", "f2-f3", "f3-f4", "f4-f5", "f5-f6")


def write_orthoxml(path, with_taxonomy=False):
    ids = {}
    species_xml = []
    for i, (sp, contigs) in enumerate(sorted(GENOMES.items())):
        genes_xml = []
        for j, name in enumerate(g for contig in contigs for g in contig):
            ids[(sp, name)] = gid = str(1000 * (i + 1) + j)
            genes_xml.append('<gene id="%s" protId="%s_%s"/>' % (gid, sp, name))
        species_xml.append('<species name="%s" NCBITaxId="%d" taxonId="%d"><database name="test" version="1"><genes>%s'
                           '</genes></database></species>' % (sp, i + 1, i + 1, "".join(genes_xml)))
    groups_xml = []
    for fam in FAMILIES:
        subgroups = "".join(
            '<orthologGroup id="%s_%s"><property name="TaxRange" value="%s"/>%s</orthologGroup>'
            % (fam, clade, clade, "".join('<geneRef id="%s"/>' % ids[(sp, fam)] for sp in clade))
            for clade in ("AB", "CD"))
        groups_xml.append('<orthologGroup id="%s"><property name="TaxRange" value="Root"/>%s</orthologGroup>'
                          % (fam, subgroups))
    path.write_text('<?xml version="1.0" encoding="utf-8"?>\n'
                    '<orthoXML xmlns="http://orthoXML.org/2011/" origin="test" originVersion="1" version="0.3">'
                    '%s%s<groups>%s</groups></orthoXML>'
                    % ("".join(species_xml), TAXONOMY if with_taxonomy else "", "".join(groups_xml)))


# taxonomy from a newick file, or from the <taxonomy> of the orthoXML without a tree file
@pytest.fixture(scope="module", params=[(False, "newick"), (True, "newick"), (False, "orthoxml"), (True, "orthoxml")],
                ids=["unoriented-newick", "oriented-newick", "unoriented-orthoxml", "oriented-orthoxml"])
def inferred(request, tmp_path_factory):
    orient, taxonomy = request.param
    tmp = tmp_path_factory.mktemp("synthetic")
    (tmp / "tree.nwk").write_text(SPECIES_TREE)
    write_orthoxml(tmp / "hogs.orthoxml", with_taxonomy=taxonomy == "orthoxml")
    tree = str(tmp / "tree.nwk") if taxonomy == "newick" else None
    gff_dir = tmp / "gff"
    gff_dir.mkdir()
    for sp, contigs in GENOMES.items():
        rows = [gene("chr%d" % (c + 1), 1000 * k + 1, 1000 * k + 500, "+", "%s_%s_gene" % (sp, name), "%s_%s" % (sp, name))
                for c, contig in enumerate(contigs) for k, name in enumerate(contig)]
        write_gff(gff_dir / ("%s.gff" % sp), *rows)

    ham = map_hogs_onto_tree(str(tmp / "hogs.orthoxml"), tree, None)
    hogxml_entries, protein_id_to_entry = get_hogxml_entries(ham)
    characterize_tree(ham.taxonomy.tree)
    init_extant_graphs(1, ham, None, str(gff_dir), hogxml_entries, protein_id_to_entry, orient)
    leaves_to_root_synteny_propagation(ham, 3, orient)
    root_to_leaves_edge_trimming(ham, orient)
    linearize_graphs(ham)
    return {node.name: node for node in ham.taxonomy.tree.traverse()}


def family(hog):
    families = {g.prot_id.split("_", 1)[1] for g in hog.get_all_descendant_genes()}
    assert len(families) == 1
    return families.pop()


def family_edges(graph):
    return {frozenset((family(a), family(b))) for a, b in graph.edges()}


def edge_between(graph, fam_a, fam_b):
    for a, b, data in graph.edges(data=True):
        if {family(a), family(b)} == {fam_a, fam_b}:
            return data
    raise KeyError((fam_a, fam_b))


def test_ancestral_gene_content(inferred):
    for level in ("AB", "CD", "Root"):
        assert sorted(family(h) for h in inferred[level].props["linear_synteny"].nodes) == FAMILIES


def test_bottom_up_contains_novel_adjacencies(inferred):
    assert edges("f2-f4", "f3-f5", "f6-f7") <= family_edges(inferred["Root"].props["bottom_up_synteny"])


def test_top_down_removes_unparsimonious_adjacencies(inferred):
    assert family_edges(inferred["Root"].props["top_down_synteny"]) == edges(*CHAIN)
    assert family_edges(inferred["CD"].props["top_down_synteny"]) == edges(*CHAIN)
    assert family_edges(inferred["AB"].props["top_down_synteny"]) == edges(*CHAIN, "f6-f7")


def test_linearized_ancestral_orders(inferred):
    assert family_edges(inferred["Root"].props["linear_synteny"]) == edges(*CHAIN)
    assert family_edges(inferred["CD"].props["linear_synteny"]) == edges(*CHAIN)
    assert family_edges(inferred["AB"].props["linear_synteny"]) == edges(*CHAIN, "f6-f7")


def test_gained_gene_is_bridged(inferred):
    # in A, x lies between f1 and f2; both A and B support f1-f2 in AB
    assert edge_between(inferred["AB"].props["bottom_up_synteny"], "f1", "f2")["weight"] == 2
    assert edge_between(inferred["Root"].props["bottom_up_synteny"], "f1", "f2")["weight"] == 4


def test_weights_count_supporting_extant_genomes(inferred):
    root = inferred["Root"].props["bottom_up_synteny"]
    assert edge_between(root, "f3", "f4")["weight"] == 4
    assert edge_between(root, "f2", "f3")["weight"] == 3
    assert sorted(c.name for c in edge_between(root, "f2", "f3")["children"]) == ["AB", "CD"]
    assert sorted(g.name for g in edge_between(root, "f2", "f3")["extant_descendants"]) == ["A", "B", "C"]


def test_linearized_graphs_have_max_degree_two(inferred):
    for level in ("AB", "CD", "Root"):
        assert max(d for _, d in inferred[level].props["linear_synteny"].degree()) <= 2
