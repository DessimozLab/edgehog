import gzip
import os
import sys
from collections import defaultdict

import pytest

from edgehog.edgehog import main

TESTS = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(os.path.dirname(TESTS), "test_data")
SNAPSHOT = os.path.join(TESTS, "snapshots", "test_data_pipeline.tsv.gz")
UPDATE_SNAPSHOT = os.environ.get("EDGEHOG_UPDATE_SNAPSHOT") == "1"


@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("edgehog_output")
    argv = ["edgehog",
            "--hogs", os.path.join(TEST_DATA, "FastOMA_HOGs.orthoxml"),
            "--species_tree", os.path.join(TEST_DATA, "species_tree.nwk"),
            "--gff_directory", os.path.join(TEST_DATA, "gff3"),
            "--date_edges",
            "--orient_edges",
            "--phylostratify",
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


def canonical_rows(out_dir):
    """All output rows per file name, in an order that does not depend on the run.

    The order of the two genes of an edge and the row order vary between runs
    (they follow set/dict iteration over pyham objects), so the genes of each
    edge are sorted and the rows of each file are sorted. The header is kept
    as the first row.
    """
    rows = {}
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.endswith(".gz"):
            with gzip.open(path, "rt") as fh:
                lines = fh.read().splitlines()
            body = []
            for line in lines[1:]:
                fields = line.split("\t")
                fields[:2] = sorted(fields[:2])
                body.append("\t".join(fields))
        else:
            with open(path) as fh:
                lines = fh.read().splitlines()
            body = lines[1:]
        rows[name] = [lines[0]] + sorted(body)
    return rows


def write_snapshot(rows):
    os.makedirs(os.path.dirname(SNAPSHOT), exist_ok=True)
    # mtime=0 so that regenerating an unchanged snapshot gives identical bytes
    with open(SNAPSHOT, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        for name, lines in rows.items():
            for line in lines:
                gz.write(("%s\t%s\n" % (name, line)).encode())


def read_snapshot():
    rows = defaultdict(list)
    with gzip.open(SNAPSHOT, "rt") as fh:
        for line in fh:
            name, row = line.rstrip("\n").split("\t", 1)
            rows[name].append(row)
    return dict(rows)


def test_output_matches_snapshot(pipeline_output):
    """Every output row on test_data matches tests/snapshots/test_data_pipeline.tsv.gz.

    After an intended change of the results, regenerate the snapshot with
    EDGEHOG_UPDATE_SNAPSHOT=1 uv run pytest tests/test_pipeline.py
    and commit it together with the change.
    """
    actual = canonical_rows(pipeline_output)
    if UPDATE_SNAPSHOT:
        write_snapshot(actual)
        pytest.skip("snapshot updated: " + SNAPSHOT)
    expected = read_snapshot()
    problems = []
    for name in sorted(set(expected) | set(actual)):
        if name not in actual:
            problems.append("%s: missing from output" % name)
        elif name not in expected:
            problems.append("%s: not in snapshot" % name)
        elif actual[name] != expected[name]:
            exp, act = set(expected[name]), set(actual[name])
            missing, extra = sorted(exp - act), sorted(act - exp)
            problems.append("%s: %d rows missing, %d extra (of %d)\n  - %s\n  + %s"
                            % (name, len(missing), len(extra), len(expected[name]),
                               "\n  - ".join(missing[:3]), "\n  + ".join(extra[:3])))
    assert not problems, "output differs from snapshot:\n" + "\n".join(problems)
