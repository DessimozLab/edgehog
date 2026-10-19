"""Tests for the species tree handling: tree characterization and the optional --species_tree."""

import argparse

import pytest
from ete4 import Tree

from edgehog.characterize_species_tree import characterize_tree
from edgehog.check_args import check_args
from edgehog.process_hogs import map_hogs_onto_tree
from test_synthetic_pipeline import write_orthoxml


def characterized(newick):
    tree = Tree(newick, parser=1)
    characterize_tree(tree)
    return {n.name: (n.props["red"], n.props["level"], n.props["n_desc_leaves"], n.props["n_internal_nodes"])
            for n in tree.traverse()}


def test_red_score_with_branch_lengths():
    nodes = characterized("((A:1,B:3)AB:2,C:1)Root;")
    # AB: d = 2, u = (sum of distances to its leaves + d * leaves) / leaves = (4 + 2 * 2) / 2 = 4
    assert nodes["AB"][0] == pytest.approx(0.5)
    assert nodes["Root"] == (0, 0, 3, 2)
    assert nodes["AB"][1:] == (1, 2, 1)
    assert nodes["A"] == (1, 2, 0, 0)


def test_missing_branch_lengths_count_as_one():
    assert characterized("((A,B)AB,(C,D)CD)Root;") == characterized("((A:1,B:1)AB:1,(C:1,D:1)CD:1)Root;")


def _args(tmp_path, species_tree):
    hogs = tmp_path / "hogs.orthoxml"
    hogs.write_text("")
    return argparse.Namespace(output_directory=str(tmp_path / "out"), species_tree=species_tree, hogs=str(hogs),
                              gff_directory=str(tmp_path), hdf5=None)


def test_species_tree_is_optional(tmp_path):
    assert check_args(_args(tmp_path, None)) == str(tmp_path / "out")


def test_missing_species_tree_file_exits(tmp_path):
    with pytest.raises(SystemExit):
        check_args(_args(tmp_path, str(tmp_path / "missing.nwk")))


def test_orthoxml_without_taxonomy_needs_species_tree(tmp_path):
    write_orthoxml(tmp_path / "hogs.orthoxml")
    with pytest.raises(SystemExit, match="pyham's error message"):
        map_hogs_onto_tree(str(tmp_path / "hogs.orthoxml"), None, None)


def test_species_tree_takes_precedence_over_orthoxml_taxonomy(tmp_path):
    write_orthoxml(tmp_path / "hogs.orthoxml", with_taxonomy=True)
    # same topology, but other internal node names than in the <taxonomy>
    (tmp_path / "tree.nwk").write_text("((A,B)X,(C,D)Y)Z;")
    ham = map_hogs_onto_tree(str(tmp_path / "hogs.orthoxml"), str(tmp_path / "tree.nwk"), None)
    assert {n.name for n in ham.taxonomy.tree.traverse()} == {"A", "B", "C", "D", "X", "Y", "Z"}
