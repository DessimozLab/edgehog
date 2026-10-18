import gzip
import os
import sys

import pytest

from edgehog.edgehog import main

TEST_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data")


@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("edgehog_output")
    argv = ["edgehog",
            "--hogs", os.path.join(TEST_DATA, "FastOMA_HOGs.orthoxml"),
            "--species_tree", os.path.join(TEST_DATA, "species_tree.nwk"),
            "--gff_directory", os.path.join(TEST_DATA, "gff3"),
            "--date_edges",
            "--orient_edges",
            "--output_directory", str(out_dir)]
    old_argv = sys.argv
    sys.argv = argv
    try:
        main()
    finally:
        sys.argv = old_argv
    return out_dir


def test_output_files_written(pipeline_output):
    files = os.listdir(pipeline_output)
    assert "genome_dict.tsv" in files
    for kind in ("bottom-up", "top-down", "linearized"):
        assert "12_{}_synteny_graph_edges.tsv.gz".format(kind) in files


def test_known_ancestral_adjacency(pipeline_output):
    expected = {"HOG_HOG:0000285_13", "HOG_HOG:0000284_13"}
    with gzip.open(os.path.join(pipeline_output, "12_linearized_synteny_graph_edges.tsv.gz"), "rt") as fh:
        matches = [line for line in fh if set(line.split("\t")[:2]) == expected]
    assert len(matches) == 1
