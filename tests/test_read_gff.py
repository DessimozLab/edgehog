"""Tests for the GFF parser (read_gff.py)."""

import pytest

from edgehog.read_gff import gff_to_dict
from gff_helpers import gene, write_gff


def gene_of(protein_dict, gene_dict, prot):
    p = protein_dict[prot]
    return gene_dict[(p["genomic_accession"], p["gene_index"])]["gene_id"]


def test_basic_gff(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 300, 400, "-", "g2", "p2"),
                    gene("chr2", 50, 90, "+", "g3", "p3"))
    gene_dict, protein_dict = gff_to_dict(gff)
    assert [v["gene_id"] for v in gene_dict.values()] == ["g1", "g2", "g3"]
    assert [k[0] for k in gene_dict] == ["chr1", "chr1", "chr2"]
    assert set(protein_dict) == {"p1", "p2", "p3"}
    assert protein_dict["p2"]["strand"] == "-"
    assert (protein_dict["p3"]["start"], protein_dict["p3"]["end"]) == (50, 90)
    for prot, gid in [("p1", "g1"), ("p2", "g2"), ("p3", "g3")]:
        assert gene_of(protein_dict, gene_dict, prot) == gid


def test_isoforms_map_to_same_gene(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1a", "p1b"),
                    gene("chr1", 300, 400, "+", "g2", "p2"))
    gene_dict, protein_dict = gff_to_dict(gff)
    assert len(gene_dict) == 2
    assert gene_of(protein_dict, gene_dict, "p1a") == gene_of(protein_dict, gene_dict, "p1b") == "g1"


def test_multi_exon_cds_gives_one_protein(tmp_path):
    rows = gene("chr1", 100, 200, "+", "g1", "p1")
    rows.append(("chr1", "CDS", 300, 400, "+", "ID=cds-g1-1;Parent=g1;protein_id=p1"))
    gff = write_gff(tmp_path / "SP.gff", rows, gene("chr1", 500, 600, "+", "g2", "p2"))
    gene_dict, protein_dict = gff_to_dict(gff)
    assert set(protein_dict) == {"p1", "p2"}
    assert protein_dict["p1"]["start"] == 100  # first CDS line is kept


def test_noncoding_and_pseudogenes_are_dropped(tmp_path):
    gff = write_gff(tmp_path / "SP.gff",
                    gene("chr1", 100, 200, "+", "g1", "p1"),
                    gene("chr1", 250, 280, "+", "nc1"),
                    gene("chr1", 300, 400, "+", "ps1", "pp1", type="pseudogene"),
                    gene("chr1", 500, 600, "+", "g2", "p2"))
    gene_dict, protein_dict = gff_to_dict(gff)
    assert [v["gene_id"] for v in gene_dict.values()] == ["g1", "g2"]
    assert set(protein_dict) == {"p1", "p2"}


def test_other_feature_types_and_comments_are_ignored(tmp_path):
    rows = [("chr1", "region", 1, 10000, "+", "ID=chr1")]
    rows += gene("chr1", 100, 200, "+", "g1", "p1")
    rows.insert(2, ("chr1", "mRNA", 100, 200, "+", "ID=rna1;Parent=g1"))
    rows.insert(3, ("chr1", "exon", 100, 200, "+", "ID=exon1;Parent=rna1"))
    gff = write_gff(tmp_path / "SP.gff", rows, [("chr1", "tRNA", 250, 300, "+", "ID=t1")],
                    gene("chr1", 300, 400, "+", "g2", "p2"),
                    header="##gff-version 3\n#!genome-build test\n")
    gene_dict, protein_dict = gff_to_dict(gff)
    assert [v["gene_id"] for v in gene_dict.values()] == ["g1", "g2"]
    assert gene_of(protein_dict, gene_dict, "p1") == "g1"


def test_feature_type_case_is_ignored(tmp_path):
    rows = [("chr1", "Gene", 100, 200, "+", "ID=g1"),
            ("chr1", "cds", 100, 200, "+", "ID=c1;protein_id=p1"),
            ("chr1", "GENE", 300, 400, "+", "ID=g2"),
            ("chr1", "Cds", 300, 400, "+", "ID=c2;protein_id=p2")]
    gene_dict, protein_dict = gff_to_dict(write_gff(tmp_path / "SP.gff", rows))
    assert set(protein_dict) == {"p1", "p2"}


def test_attribute_parsing(tmp_path):
    rows = [("chr1", "gene", 100, 200, "+", "ID=g1;Note=a=b;"),
            ("chr1", "CDS", 100, 200, "+", "ID=c1;Parent=g1;product=x=y;protein_id=p1;")]
    gene_dict, protein_dict = gff_to_dict(write_gff(tmp_path / "SP.gff", rows))
    assert list(protein_dict) == ["p1"]
    assert [v["gene_id"] for v in gene_dict.values()] == ["g1"]


def test_gene_id_preferred_over_id(tmp_path):
    rows = [("chr1", "gene", 100, 200, "+", "ID=gene-1;gene_id=ENSG1"),
            ("chr1", "CDS", 100, 200, "+", "ID=c1;Parent=gene-1;protein_id=p1")]
    gene_dict, _ = gff_to_dict(write_gff(tmp_path / "SP.gff", rows))
    assert [v["gene_id"] for v in gene_dict.values()] == ["ENSG1"]


def test_without_gene_or_cds_lines_exits(tmp_path):
    gff = write_gff(tmp_path / "SP.gff", [("chr1", "region", 1, 100, "+", "ID=chr1")])
    with pytest.raises(SystemExit):
        gff_to_dict(gff)


def test_malformed_gff_exits(tmp_path):
    path = tmp_path / "SP.gff"
    path.write_text("chr1\ttest\tgene\tstart\tend\t.\t+\t.\tID=g1\n")
    with pytest.raises(SystemExit):
        gff_to_dict(str(path))


# Issue #8 (found while writing these tests): CDS-only GFFs without
# protein_id, e.g. Prodigal output, are meant to be supported (read_gff.py:60,
# 119-128), but the derived protein_id is only added to df_proteins, so
# df['protein_id'] raises a KeyError.

@pytest.mark.xfail(strict=True, raises=KeyError, reason="#8: CDS-only (Prodigal) GFFs without protein_id crash")
def test_prodigal_cds_only_gff(tmp_path):
    rows = [("contigA", "CDS", 1, 300, "+", "ID=1_1;partial=00"),
            ("contigA", "CDS", 400, 900, "-", "ID=1_2;partial=00")]
    gene_dict, protein_dict = gff_to_dict(write_gff(tmp_path / "SP.gff", rows))
    assert set(protein_dict) == {"contigA_1", "contigA_2"}
    assert len(gene_dict) == 2


# Issue #9 (found while writing these tests): GFF3 files may end with an
# embedded ##FASTA section; its sequence lines are parsed as features.

@pytest.mark.xfail(strict=True, raises=SystemExit, reason="#9: GFFs with an embedded ##FASTA section are rejected")
def test_embedded_fasta_section(tmp_path):
    path = write_gff(tmp_path / "SP.gff", gene("chr1", 100, 200, "+", "g1", "p1"))
    with open(path, "a") as fh:
        fh.write("##FASTA\n>chr1\nACGTACGT\n")
    _, protein_dict = gff_to_dict(path)
    assert set(protein_dict) == {"p1"}
