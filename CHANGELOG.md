# Changelog

## 0.2.0 (unreleased)

This release fixes several deviations from the algorithm described in the
Supplementary Information of the edgeHOG paper
([doi:10.1038/s41559-025-02818-0](https://doi.org/10.1038/s41559-025-02818-0))
and bugs in reading the input (tracked in
[#43](https://github.com/DessimozLab/edgehog/issues/43)). **The inferred
ancestral synteny graphs change for existing datasets**, mostly because of
the gene order of unsorted GFF files.

### Changes in results
- Extant genes are ordered by contig and genomic coordinates before adjacent
  genes are linked. Before, the order of the GFF file was used, which links
  wrong genes if the GFF is not sorted. On `test_data`, about 27% of the
  linearized ancestral edges change.
  ([#40](https://github.com/DessimozLab/edgehog/issues/40),
  [#46](https://github.com/DessimozLab/edgehog/pull/46))
- The orientation of an extant edge uses the strand of the previous linked
  gene, not of a gene skipped because it is missing from the orthoXML.
  ([#40](https://github.com/DessimozLab/edgehog/issues/40),
  [#46](https://github.com/DessimozLab/edgehog/pull/46))
- Bottom-up phase: if the edge reconnecting the neighbors of a gene without
  parent HOG already exists, the minimal weight of the two removed edges is
  added to it, as described in the paper. Before, it was left unchanged, and
  the two removed edges could be pruned in the top-down phase. For gaps of
  three or more genes, all real edges are now mapped to their reconnected
  edge, independent of the order in which the genes are removed.
  ([#38](https://github.com/DessimozLab/edgehog/issues/38),
  [#51](https://github.com/DessimozLab/edgehog/pull/51))
- Top-down phase with `--orient_edges`: orientation contexts not supported
  by two children or by the parent are reset to 0, as described in the paper;
  if none is supported, they stay unchanged. Before, the bottom-up contexts
  were kept. `nb_internal_nodes_from_ancestor_with_updated_weight` no longer
  depends on `--orient_edges`.
  ([#36](https://github.com/DessimozLab/edgehog/issues/36),
  [#49](https://github.com/DessimozLab/edgehog/pull/49))
- Linearization: removed edges are reintroduced by descending weight, then by
  descending weight of the resulting connected component (the second
  criterion was inverted), and no longer when they would close a cycle.
  Ties are broken in a deterministic order.
  ([#37](https://github.com/DessimozLab/edgehog/issues/37),
  [#50](https://github.com/DessimozLab/edgehog/pull/50))
- With pyham 2, ancestral genes follow the HOGs of the orthoXML. pyham 1
  added spurious ancestral genes without edges for some families with
  duplications. Some HOG labels in the output change, as pyham 2 gives HOGs
  of paralogous branches their own ids.
  ([#48](https://github.com/DessimozLab/edgehog/pull/48))

### Added
- `--species_tree` is optional if the orthoXML contains a `<taxonomy>`
  section (e.g. orthoXML exported from the OMA browser). If both are given,
  the species tree is used.
  ([#48](https://github.com/DessimozLab/edgehog/pull/48))
- Support for Python 3.13 and 3.14. The `oma` extra (needed for `--hdf5`)
  is only available for Python < 3.13, because pyoma requires numpy < 2.
- CDS-only GFF files without `protein_id` (e.g. Prodigal output) and GFF3
  files with a `##FASTA` section can be read. For Prodigal IDs, the protein
  ID is `<contig>_<orf number>`, as in the headers of Prodigal's protein
  fasta. ([#41](https://github.com/DessimozLab/edgehog/issues/41),
  [#45](https://github.com/DessimozLab/edgehog/pull/45))

### Changed
- Depends on pyham 2 and ete4 instead of pyham 1 and ete3. ete4 is only
  distributed as source, so installing edgehog needs a C compiler.
  ([#48](https://github.com/DessimozLab/edgehog/pull/48))
- GFF lookup: each GFF file is attributed to the longest genome name it
  starts with, and `<name>.gff` / `<name>.gff3` is preferred. Only files
  ending in `.gff` or `.gff3`, optionally compressed (`.gz`, `.bz2`, `.xz`,
  `.zip`), are considered. If several files remain for a genome, edgehog exits
  with an error listing them; before, it used the first one in directory
  order. ([#42](https://github.com/DessimozLab/edgehog/issues/42),
  [#44](https://github.com/DessimozLab/edgehog/pull/44))
- For CDS lines with an `ID` but no `protein_id`, the ID is used as protein
  ID (unless it is in Prodigal's format). Before, it was truncated to
  `<contig>_<second part of the ID>`, which failed for IDs without `_`.
  ([#45](https://github.com/DessimozLab/edgehog/pull/45))
- In the output TSV files, the two genes of an edge can be written in either
  order, and `weight` is written as an integer when a graph has no singletons.

### Fixed
- A missing GFF file exits with an error message instead of raising
  `IndexError`. ([#44](https://github.com/DessimozLab/edgehog/pull/44))
- A genome whose name is a prefix of another genome's name (`HUMAN` and
  `HUMAN2`) could read the other genome's GFF.
  ([#44](https://github.com/DessimozLab/edgehog/pull/44))

### Development
- Build and dependency management switched from poetry to uv.
  ([#34](https://github.com/DessimozLab/edgehog/pull/34))
- Unit test suite for GFF parsing and every inference phase, including the
  example of Fig. 2 of the paper.
  ([#35](https://github.com/DessimozLab/edgehog/pull/35))
