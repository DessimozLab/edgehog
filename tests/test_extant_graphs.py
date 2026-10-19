"""Tests for building extant synteny graphs from GFFs (init_extant_synteny_graphs.py)."""

import os
import pytest

from edgehog import init_extant_synteny_graphs as extant
from gff_helpers import gene, write_gff
from synteny_helpers import Gene, edge_set, edges


class Genome:
    def __init__(self, name):
        self.name = name


class StubHam:
    """Maps orthoXML entry ids to genes; a gene's HOG is its parent."""

    def __init__(self, entry_to_gene):
        self.entry_to_gene = entry_to_gene

    def get_gene_by_id(self, entry):
        return self.entry_to_gene[entry]

    def get_hog_by_gene(self, gene):
        return gene.parent


def build_graph(gff, in_orthoxml, orient_edges=True, prot_id_format="{}"):
    """Run gff_to_graph for genome SP; ``in_orthoxml`` lists the proteins in the orthoXML."""
    genome = Genome("SP")
    genes = {p: Gene(p, Gene("HOG_" + p)) for p in in_orthoxml}
    ham = StubHam({"e_" + p: g for p, g in genes.items()})
    prot_to_entry = {genome: {prot_id_format.format(p): "e_" + p for p in in_orthoxml}}
    graph, _ = extant.gff_to_graph(ham, genome, gff, {genome: []}, prot_to_entry, orient_edges)
    return graph


def context(graph, a, b):
    by_name = {str(n): n for n in graph}
    data = graph[by_name[a]][by_name[b]]
    return {k for k in ("unidirectional", "convergent", "divergent") if data[k]}


def test_chain_without_edges_across_contigs(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 300, 400, "+", "g2", "p2"),
                    gene("chr1", 500, 600, "+", "g3", "p3"),
                    gene("chr2", 100, 200, "+", "g4", "p4"))
    graph = build_graph(gff, ["p1", "p2", "p3", "p4"])
    assert edge_set(graph) == edges("p1-p2", "p2-p3")
    assert "p4" in {str(n) for n in graph}
    assert all(d["weight"] == 1 for _, _, d in graph.edges(data=True))
    assert {graph.nodes[n]["contig"] for n in graph} == {"chr1", "chr2"}


def test_genes_missing_from_orthoxml_are_skipped(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 50, 90, "+", "g0", "p0"),
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 300, 400, "+", "gx", "px"),
                    gene("chr1", 500, 600, "+", "g2", "p2"))
    graph = build_graph(gff, ["p1", "p2"])
    assert edge_set(graph) == edges("p1-p2")


def test_transcriptional_orientation(tmp_path):
    strands = ["+", "+", "-", "+", "-", "-"]
    gff = write_gff(tmp_path / "SP.gff",
                    *[gene("chr1", 100 * i + 1, 100 * i + 50, s, "g%d" % i, "p%d" % i) for i, s in enumerate(strands)])
    graph = build_graph(gff, ["p%d" % i for i in range(len(strands))])
    assert context(graph, "p0", "p1") == {"unidirectional"}  # -> ->
    assert context(graph, "p1", "p2") == {"convergent"}      # -> <-
    assert context(graph, "p2", "p3") == {"divergent"}       # <- ->
    assert context(graph, "p4", "p5") == {"unidirectional"}  # <- <-


def test_no_orientation_without_option(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 300, 400, "+", "g2", "p2"))
    graph = build_graph(gff, ["p1", "p2"], orient_edges=False)
    assert context(graph, "p1", "p2") == set()


@pytest.mark.parametrize("fmt", ["{} some description", "{},extra"])
def test_orthoxml_protein_ids_with_suffix(tmp_path, fmt):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 300, 400, "+", "g2", "p2"))
    graph = build_graph(gff, ["p1", "p2"], prot_id_format=fmt)
    assert edge_set(graph) == edges("p1-p2")


def test_representative_hog_is_majority_among_isoforms():
    hog_a, hog_b = Gene("A"), Gene("B")
    ham = StubHam({"e1": Gene("i1", hog_a), "e2": Gene("i2", hog_b), "e3": Gene("i3", hog_b)})
    assert extant.get_representative_hog_for_gene(ham, ["e1", "e2", "e3"]).parent is hog_b


def test_find_gff_files(tmp_path):
    for name in ["HUMAN.gff", "MOUSE.v2.gff3", "README.txt"]:
        (tmp_path / name).write_text("")
    genomes = [Genome("HUMAN"), Genome("MOUSE")]
    found = extant.get_gffs_of_extant_genomes(str(tmp_path), genomes)
    assert {g.name: os.path.basename(p) for g, p in found.items()} == {"HUMAN": "HUMAN.gff", "MOUSE": "MOUSE.v2.gff3"}


# Issue #6: the gene order is taken from the GFF file order; nothing sorts the
# features by contig and start coordinate (read_gff.py:100-117,
# init_extant_synteny_graphs.py:88).

@pytest.mark.xfail(strict=True, reason="#6: unsorted GFFs give wrong adjacencies")
def test_unsorted_gff(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 500, 600, "+", "g3", "p3"),
                    gene("chr1", 300, 400, "+", "g2", "p2"))
    assert edge_set(build_graph(gff, ["p1", "p2", "p3"])) == edges("p1-p2", "p2-p3")


@pytest.mark.xfail(strict=True, reason="#6: contigs interleaved in the GFF lose adjacencies")
def test_interleaved_contigs(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr2", 100, 200, "+", "h1", "q1"),
                    gene("chr1", 300, 400, "+", "g2", "p2"))
    assert edge_set(build_graph(gff, ["p1", "p2", "q1"])) == edges("p1-p2")


# Issue #7: the orientation uses the strand of gene_keys[i-1], which is not the
# previous linked gene if a gene missing from the orthoXML lies in between
# (init_extant_synteny_graphs.py:111).

@pytest.mark.xfail(strict=True, reason="#7: orientation uses the strand of a skipped gene")
def test_orientation_ignores_skipped_gene(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 300, 400, "-", "gx", "px"),
                    gene("chr1", 500, 600, "+", "g2", "p2"))
    assert context(build_graph(gff, ["p1", "p2"]), "p1", "p2") == {"unidirectional"}


# Found while writing these tests: a genome name that is a prefix of another
# one can match the other genome's GFF (regex '^NAME.*\.gff', first match in
# os.listdir order), and a missing GFF raises IndexError instead of the
# intended error message (init_extant_synteny_graphs.py:18-22).

@pytest.mark.xfail(strict=True, reason="genome name prefix matches another genome's GFF")
def test_find_gff_files_prefix_ambiguity(tmp_path, monkeypatch):
    monkeypatch.setattr(extant.os, "listdir", lambda d: ["HUMAN2.gff", "HUMAN.gff"])
    found = extant.get_gffs_of_extant_genomes(str(tmp_path), [Genome("HUMAN")])
    assert os.path.basename(next(iter(found.values()))) == "HUMAN.gff"


@pytest.mark.xfail(strict=True, raises=IndexError, reason="missing GFF raises IndexError instead of exiting with a message")
def test_find_gff_files_missing(tmp_path):
    with pytest.raises(SystemExit):
        extant.get_gffs_of_extant_genomes(str(tmp_path), [Genome("HUMAN")])
