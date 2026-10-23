# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

edgeHOG infers ancestral gene order (synteny graphs) at every internal node of a species tree, using HOGs (Hierarchical Orthologous Groups) from an orthoXML file plus extant genome annotations (GFF directory or an OMA HDF5 database). Python >=3.9,<3.15, managed with uv (build backend `uv_build`, flat layout via `module-root = ""`); CLI entry point `edgehog = edgehog.edgehog:main`.

## Commands

```bash
uv sync --extra oma             # dev install incl. the default `dev` group; the `oma` extra (pyoma + pytables) is needed for --hdf5 input/output
uv run edgehog --help
uv lock                         # after changing dependencies; CI uses `uv sync --locked`, so uv.lock must be committed and current

# Reference end-to-end run (the same run tests/test_pipeline.py does)
uv run edgehog --hog test_data/FastOMA_HOGs.orthoxml \
    --species_tree test_data/species_tree.nwk \
    --gff_directory test_data/gff3/ \
    --date_edges --orient_edges \
    --output_directory test-results

# Check a known ancestral edge by hand: exactly one matching line must exist
gunzip -c test-results/12_linearized_synteny_graph_edges.tsv.gz \
  | awk '/(HOG_HOG:0000285_13\tHOG_HOG:0000284_13)|(HOG_HOG:0000284_13\tHOG_HOG:0000285_13)/{c++} END{exit c==1?0:1}'
```

```bash
uv run pytest                   # ~40s; tests/test_pipeline.py runs the test_data pipeline above in-process
uv run pytest tests/test_pipeline.py::test_known_ancestral_adjacency
uv run --python 3.9 pytest      # check the oldest supported Python
```

`tests/test_pipeline.py` runs the test_data pipeline above once per module and checks the same edge. Unit tests per phase (`test_read_gff.py`, `test_extant_graphs.py`, `test_bottom_up.py`, `test_top_down.py`, `test_linearize.py`) use stubs from `tests/synteny_helpers.py` (genes with `.parent`, an ete4 tree with fake pyham genomes patched in via `patch_pyham`) instead of real orthoXML; `test_synthetic_pipeline.py` runs the phases with real pyham on a hand-made 4-species dataset with known ancestral orders. Known deviations from the paper's Supplementary Information and other bugs are documented as `xfail(strict=True)` tests, so a fix makes them fail until the marker is removed. CI (`.github/workflows/run-tests.yml`, Python 3.9–3.14) runs them after `uv run edgehog --version`, which checks the console-script entry point. An `experimental` matrix entry (currently 3.15) tests the next Python version with `continue-on-error`: it widens `requires-python` in its checkout to the matrix version plus one minor and runs `uv sync` without `--locked`; once its dependencies have wheels, move it into the main version list and bump `requires-python`. The `oma` extra is limited to Python <3.13 because pyoma requires numpy<2, which has no wheels for 3.13+, so `--hdf5` only works there. `test_data/oma-hdf5-input/` contains a small OMA HDF5 dataset (`OmaServer.h5`, `edgehog-hogs.orthoxml.gz`, `simplified_speciestree.nwk`) for exercising the `--hdf5` code path. `utils/extract_age_timetree.py` is a standalone post-processing script (dates extant edges in million years using a TimeTree tree).

## Architecture

The pipeline is orchestrated linearly in `edgehog/edgehog.py:main()`. All state lives on a single **pyham `Ham` object**: each stage stores its networkx graphs as ete4 props on the nodes of `ham.taxonomy.tree` (pyham 2 / ete4: write with `tree_node.add_prop(name, value)`, read values edgehog always sets with `tree_node.props[name]`, and the genome with `tree_node.props.get('genome')`, which is None for nodes without genome; attribute access and setting arbitrary attributes on ete4 nodes do not work, `is_leaf` is a property, and missing branch lengths are `None`, which `characterize_species_tree.branch_length` counts as 1.0). Graph nodes are pyham gene/HOG objects (extant genes at leaves, ancestral HOGs at internal nodes), so `gene.parent` maps a node to its HOG in the parent level.

Stages and the tree-node features they produce:

1. `process_hogs.map_hogs_onto_tree` — builds the `Ham` from tree + orthoXML. Whether pyham uses internal node names depends on whether the newick root has a non-numeric name; with `--hdf5`, OMA species resolve mode is used. Without `--species_tree`, pyham builds the tree from the `<taxonomy>` section of the orthoXML (if both exist, the newick wins).
2. `characterize_species_tree.characterize_tree` — computes per-node stats (RED score, level from root) used in output/genome_dict.tsv.
3. Extant graphs → `bottom_up_synteny` and `contiguity_dict` on leaves:
   - GFF path: `init_extant_synteny_graphs.init_extant_graphs` with `read_gff.py`; orthoXML `protId` (or its prefix before a space) must match the GFF CDS `protein_id`, and GFF files must be named after tree leaves.
   - HDF5 path: `hdf5.init_extant_graphs_from_hdf5` reads OMA main isoforms sorted by chromosome/locus. Genes absent from the orthoXML are skipped (including the first one). Also sets `genome_code` and `taxid`.
   - Adjacent genes are only linked when on the same contig. Edge attrs: `weight`, `unidirectional`/`convergent`/`divergent` (orientation counts, only filled with `--orient_edges`), `children`, `extant_descendants`.
4. `infer_ancestral_synteny_graphs.leaves_to_root_synteny_propagation` (postorder) — projects each child's graph onto parent HOGs, summing weights; `build_transitive_graph` bridges over up to `--max_gaps` child genes that have no parent HOG. Sets `bottom_up_synteny` on ancestors and `transitive_synteny` / `real_edge_to_transitive_edge` on children.
5. `root_to_leaves_edge_trimming` (preorder) — removes edges not supported parsimoniously (seen in only one child lineage) → `top_down_synteny`.
6. `add_ons.date_edges` (optional, `--date_edges`), then `linearize_graphs` → `linear_synteny`, then `add_ons.phylostratify` (optional) → `stratigraphy`.
7. `write_output.write_output` writes per-node gzipped TSVs (`<id>_bottom-up|top-down|linearized|extant_synteny_graph_edges.tsv.gz`) plus `genome_dict.tsv` mapping numeric prefixes to tree nodes; `write_as_hdf5` (`--out-format HDF5`) uses `hdf5.HDF5Writer`, which requires the OMA HDF5 DB to map level names to NCBI taxids and HOG IDs to row numbers.

`hdf5.py` imports pyoma/tables at module level and exits if missing, so import it lazily (as `init_extant_synteny_graphs.py` does) to keep the GFF path working without the `oma` extra.

Version is read from installed package metadata (`edgehog/__init__.py`), so bump it in `pyproject.toml`.
