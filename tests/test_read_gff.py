"""Tests for the GFF parser (read_gff.py)."""

import bz2
import gzip
import lzma
import zipfile

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


# Issue #41: CDS-only GFFs without protein_id (e.g. Prodigal output) and GFF3
# files with an embedded ##FASTA section could not be read.

PRODIGAL_GFF = """##gff-version  3
# Sequence Data: seqnum=1;seqlen=5000;seqhdr="contigA some description"
# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=50.00;transl_table=11;uses_sd=1
contigA\tProdigal_v2.6.3\tCDS\t1\t300\t45.2\t+\t0\tID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.52;conf=100.00;score=45.21;
contigA\tProdigal_v2.6.3\tCDS\t400\t900\t60.1\t-\t0\tID=1_2;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.48;conf=99.99;score=60.12;
# Sequence Data: seqnum=2;seqlen=3000;seqhdr="contigB"
# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=50.00;transl_table=11;uses_sd=1
contigB\tProdigal_v2.6.3\tCDS\t10\t500\t30.0\t+\t0\tID=2_1;partial=00;start_type=GTG;rbs_motif=None;rbs_spacer=None;gc_cont=0.50;conf=99.00;score=30.01;
"""


def test_prodigal_cds_only_gff(tmp_path):
    path = tmp_path / "SP.gff"
    path.write_text(PRODIGAL_GFF)
    gene_dict, protein_dict = gff_to_dict(str(path))
    # protein ids as in the headers of prodigal's protein fasta
    assert list(protein_dict) == ["contigA_1", "contigA_2", "contigB_1"]
    assert [(k[0], v["gene_id"], v["strand"]) for k, v in gene_dict.items()] == [
        ("contigA", "contigA_1", "+"), ("contigA", "contigA_2", "-"), ("contigB", "contigB_1", "+")]
    for prot in protein_dict:
        assert gene_of(protein_dict, gene_dict, prot) == prot


def test_cds_only_gff_with_other_ids(tmp_path):
    rows = [("chr1", "CDS", 1, 300, "+", "ID=cds-XP_001"),
            ("chr1", "CDS", 400, 900, "+", "ID=geneA"),
            ("chr1", "CDS", 1000, 1200, "+", "Name=no_id")]
    _, protein_dict = gff_to_dict(write_gff(tmp_path / "SP.gff", rows))
    assert list(protein_dict) == ["cds-XP_001", "geneA"]


def test_cds_without_protein_id_in_gff_with_genes(tmp_path):
    rows = [("chr1", "gene", 1, 300, "+", "ID=g1"),
            ("chr1", "CDS", 1, 300, "+", "ID=cds1;Parent=g1"),
            ("chr1", "gene", 400, 900, "+", "ID=g2"),
            ("chr1", "CDS", 400, 900, "+", "ID=cds2;Parent=g2")]
    gene_dict, protein_dict = gff_to_dict(write_gff(tmp_path / "SP.gff", rows))
    assert list(protein_dict) == ["cds1", "cds2"]
    assert gene_of(protein_dict, gene_dict, "cds2") == "g2"


def test_cds_only_gff_without_ids_exits(tmp_path):
    rows = [("chr1", "CDS", 1, 300, "+", "Name=a"), ("chr1", "CDS", 400, 900, "+", "Name=b")]
    with pytest.raises(SystemExit, match="Neither the tag 'protein_id' nor 'ID'"):
        gff_to_dict(write_gff(tmp_path / "SP.gff", rows))


def test_embedded_fasta_section(tmp_path):
    path = write_gff(tmp_path / "SP.gff", gene("chr1", 100, 200, "+", "g1", "p1"))
    with open(path, "a") as fh:
        fh.write("##FASTA\n>chr1\nACGTACGT\n>p1\nMKV\n")
    _, protein_dict = gff_to_dict(path)
    assert set(protein_dict) == {"p1"}


def _compress(path, suffix):
    data = path.read_bytes()
    target = path.with_name(path.name + suffix)
    if suffix == ".zip":
        with zipfile.ZipFile(target, "w") as archive:
            archive.writestr(path.name, data)
    else:
        opener = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}[suffix]
        with opener(target, "wb") as fh:
            fh.write(data)
    return str(target)


@pytest.mark.parametrize("suffix", [".gz", ".bz2", ".xz", ".zip"])
def test_compressed_gff(tmp_path, suffix):
    path = write_gff(tmp_path / "SP.gff",
                     gene("chr1", 100, 200, "+", "g1", "p1"),
                     gene("chr1", 300, 400, "-", "g2", "p2"))
    with open(path, "a") as fh:
        fh.write("##FASTA\n>chr1\nACGTACGT\n")
    gene_dict, protein_dict = gff_to_dict(_compress(tmp_path / "SP.gff", suffix))
    assert list(protein_dict) == ["p1", "p2"]
    assert gene_of(protein_dict, gene_dict, "p2") == "g2"


def test_zip_with_several_files_exits(tmp_path):
    target = tmp_path / "SP.gff.zip"
    with zipfile.ZipFile(target, "w") as archive:
        archive.writestr("a.gff", "")
        archive.writestr("b.gff", "")
    with pytest.raises(SystemExit):
        gff_to_dict(str(target))
