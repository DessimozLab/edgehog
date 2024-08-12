# edgeHOG

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation)
- [Demo](#demo)
- [Results](#results)
- [License](./LICENSE.txt)
- [Issues](https://github.com/dessimozlab/edgehog/issues)


## Overview

```edgeHOG``` is a tool to infer the gene order of each ancestor in a species phylogeny. As such, ```edgeHOG``` enables both to explore ancestral microsyntenies (local scale) and to reconstruct ancestral chromosomes (global scale). 

```edgeHOG``` relies on objects called HOGs (Hierarchical Groups of Orthologs) to model gene lineages and ancestral gene content. Basically, genes that belong to the same HOG across extant genomes are inferred to have descended from the same common ancestral gene in the common ancestor of these genomes. Accordingly, adjacencies between extant genes can be converted to edges between HOGs, which enables parsimonious ancestral gene order inferences.  

## System Requirements


### Hardware Requirements

The `edgeHOG` package requires only a standard computer with enough RAM. The amount of RAM depends a lot on the size of the dataset. For big datasets (thousands of genomes), more than 100GB of RAM are needed.

### Software Requirements

#### OS Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on an Ubuntu 22.04 and CentOS 7 environment.

The package itself should be compatible with Windows, Mac and Linux operating systems.

Edgehog is written in purge Python, so a working python installation is needed before installing edgehog.

#### Installing Python on Ubuntu 22.04

Python can be installed directly from its `apt` system using `apt install python3`

## Installation

### From PyPi using pip
`edgeHOG`  can be installed directly from pypi using pip. The command is the following:

```bash
pip install edgehog
```

The dependencies are version pinned and will automatically be installed as well. 

### From sources
```edgeHOG``` was built and tested with python 3.9 and higher. To set up ```edgeHOG``` on your local machine, please follow the instructions below. 

```bash
pip install poetry  # poetry is used as build and dependency resolving system.

git clone https://github.com/dessimozlab/edgehog.git
cd edgehog
poetry install
```

## Usage

```
usage: edgehog [-h] [--version] [--output_directory OUTPUT_DIRECTORY] 
               --species_tree SPECIES_TREE  --hogs HOGS 
               [--gff_directory GFF_DIRECTORY] [--hdf5 HDF5]  [--date_edges] 
               [--phylostratify] [--max_gaps MAX_GAPS] [--include_extant_genes] 
               [--out-format {TSV,HDF5}]

edgehog is a software tool that infers an ancestral synteny graph at each
internal node of an input species phylogenetic tree

optional arguments:
  -h, --help            show this help message and exit
  --version             print version number and exit
  --output_directory OUTPUT_DIRECTORY
                        path to output directory (default is ./edgehog_output)
  --species_tree SPECIES_TREE
                        path to species/genomes phylogenetic tree (newick format)
  --hogs HOGS           path to the HierarchicalGroups.orthoxml file in which HOGs are stored
  --gff_directory GFF_DIRECTORY
                        path to directory with the gffs of extant genomes (each gff file must be named according to the name of an extant genome / leaf on the species tree)
  --hdf5 HDF5           path to the hdf5 file (alternative to gff_directory to run edgeHOG on the entire OMA database)
  --date_edges          whether the age of edges in extant species should be predicted
  --phylostratify       whether the number of edge retention, gain and loss should be analyzed for each node of the species tree
  --max_gaps MAX_GAPS   max_gaps can be seen as the theoritical maximal number of consecutive novel genes that can emerge between two older genes (default = 3), e.g. if max_gaps = 2: the probabilistic A-b-c-D-E-f-g-h-I-J graph will be turn into A-D-E ; I-J in the ancestorwhile if max_gaps = 3: the probabilistic A-b-c-D-E-f-g-h-I-J graph will be
                        turn into A-D-E-I-J in the ancestor
  --include_extant_genes
                        include extant genes in output file for ancestral reconstructions.
  --out-format {TSV,HDF5}
                        define output format. Can be TSV (tab seperated files) or HDF5 (compatible for integration into oma hdf5)
```

## Input data

Three types of input data are needed for ```edgeHOG``` to run:
* a phylogenetic tree of species/genomes of interest (in newick format)
* the annotation of each of these genomes (e.g. in the form of a directory of gff files)
* an HierarchicalGroups.orthoxml file corresponding to the extant genomes

Since these input data are intersected, they must comply with the following requirements:
* the prefix of a gff filename must correspond to a species/genome identifier in the phylogenetic tree
* all genome identifiers in the phylogenetic tree must correspond to a genome entry in the HierarchicalGroups.orthoxml file
* the ```protId``` or the prefix of the ```protId``` followed by the ```' '``` character of each entry of a given input genome in the HierarchicalGroups.orthoxml file must match the ```protein_id``` of a CDS in the gff file of this genome

### Species tree

A phylogenetic tree of the input genomes/species must be provided in the newick format. If internal nodes are not named, they will be named based on the concatenation of the names of their descendant leaves

* The species phylogenetic tree used in the OMA database can be downloaded [here](https://omabrowser.org/All/speciestree.nwk)
* The high-quality GTDB archaeal and bacterial species trees, along with metadata can be found [here](https://data.gtdb.ecogenomic.org/releases/latest/)
* To use the tree or a subtree of the NCBI taxonomy database, the ete3 python package has some [useful build-in functions](http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html)

To prune a tree in order to obtain only the phylogeny of your genomes of interest, please refer to the corresponding [ete3 tutorial](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#pruning-trees)

If you don't have a species tree available for your genomes, you can follow [this tutorial](https://github.com/DessimozLab/f1000_PhylogeneticTree) on how to use OMA Browser and OMA standalone for species tree inference.


### HierarchicalGroups.orthoxml
 
An HierarchicalGroups.orthoxml file of HOGs defined based on the proteomes and the species tree of input genomes is required for genomic comparisons. HOGs xml files can be retrieved from several leading orthology databases such as [OrthoDB](https://www.orthodb.org/), [EggNOG](http://eggnog5.embl.de), [HieranoiDB](https://hieranoidb.sbc.su.se/) or [OMA](https://omabrowser.org/oma/home/).

If you don't have a HierarchicalGroups.orthoxml for your genomes, HOGs can be inferred from your input dataset using [OMA_standalone](https://omabrowser.org/standalone/). 

## Demo

### Small test dataset

We provide a small testdata set in the subdirectory `test_data`. `edgeHOG` can be run on this dataset with the following command:

```bash
edgehog --hog test_data/FastOMA_HOGs.orthoxml \
                --species_tree test_data/species_tree.nwk \
                --gff_directory test_data/gff3/ \
                --date_edges \
                --output_directory test-results
```

See the [test-data specific README](test_data/README.md) for more details how the dataset was assembled. The [Result section](#results) will discuss what the result files contain and how they can be interpreted.  

### Large dataset (complete OMA database with thousands of genomes)
`edgehog` can be run on the complete public OMA database using the data available on https://omabrowser.org/oma/current/. For that,
one can download the HOGs ([oma-hogs.orthoXML.gz file](https://omabrowser.org/All/oma-hogs.orthoXML.gz)), the [species tree](https://omabrowser.org/All/speciestree.nwk) and 
the [OMA HDF5 database](https://omabrowser.org/All/OmaServer.h5). 

Note that this dataset is very large (>200 GB). It can be run with the following command:

```bash
wget https://omabrowser.org/All/oma-hogs.orthoXML.gz
wget https://omabrowser.org/All/speciestree.nwk
wget https://omabrowser.org/All/OmaServer.h5

gunzip oma-hogs.orthoXML.gz
edgehog --hogs oma-hogs.orthoXML --hdf5 OmaServer.h5 --species_tree speciestree.nwk --date_edges --output_directory ./edghog_results
```



## Results

edgehog produces a number of result files in the specified output directory (e.g. `./edgehog_output`). Unless the `--out-format` is 
specified to be hdf5, the result files are all TSV files:

```bash
$> ls edgehog_results
0_bottom-up_synteny_graph_edges.tsv.gz   4_extant_synteny_graph_edges.tsv.gz
0_linearized_synteny_graph_edges.tsv.gz  5_extant_synteny_graph_edges.tsv.gz
0_top-down_synteny_graph_edges.tsv.gz    6_bottom-up_synteny_graph_edges.tsv.gz
1_bottom-up_synteny_graph_edges.tsv.gz   6_linearized_synteny_graph_edges.tsv.gz
1_linearized_synteny_graph_edges.tsv.gz  6_top-down_synteny_graph_edges.tsv.gz
1_top-down_synteny_graph_edges.tsv.gz    7_extant_synteny_graph_edges.tsv.gz
2_bottom-up_synteny_graph_edges.tsv.gz   8_extant_synteny_graph_edges.tsv.gz
2_linearized_synteny_graph_edges.tsv.gz  9_bottom-up_synteny_graph_edges.tsv.gz
2_top-down_synteny_graph_edges.tsv.gz    9_linearized_synteny_graph_edges.tsv.gz
3_bottom-up_synteny_graph_edges.tsv.gz   9_top-down_synteny_graph_edges.tsv.gz
3_linearized_synteny_graph_edges.tsv.gz  genome_dict.tsv
3_top-down_synteny_graph_edges.tsv.gz
```

The genome_dict.tsv file will provide a mapping from the species tree nodes to the prefix of the result files:

```TSV
genome_id       nb_descendant_leaves    level_from_root RED_score       name
0       85      0       0.00    Viridiplantae
1       7       1       0.26    Chlorophyta
2       4       2       0.51    Mamiellales
3       2       3       0.75    Ostreococcus
4       0       4       1.00    Ostreococcus tauri
5       0       4       1.00    Ostreococcus lucimarinus (strain CCE9901)
6       2       3       0.75    Micromonas
7       0       4       1.00    Micromonas commoda (strain RCC299 / NOUM17 / CCMP2709)
8       0       4       1.00    Micromonas pusilla (strain CCMP1545)
9       3       2       0.54    core chlorophytes
```

Files starting with `0_` will therefor be describing the synteny at the level of Viridiplantae (the root node of this dataset). We can see that this taxonomic level contains 85 species in total.

For internal taxonomic levels (ancestral nodes), edgehog produces three TSV files each, one for the bottom up phase 
where extant adjacencies are propagated and collected (`0_bottom-up_synteny_graph_edges.tsv.gz`), one for the 
top-down phase edges where non-parsimonious edges are removed (`0_top-down_synteny_graph_edges.tsv.gz`) and one which 
contains a linearized form (subset of top-down) that corresponds to our proposed ancestral order 
(`0_linearized_synteny_graph_edges.tsv.gz`).

Those files contain the following columns:
 - gene1: extant/ancestral gene-id of the first gene.
 - gene2: extant/ancestral gene-id of the second gene.    
 - weight: number of extant edges supporting this adjacency
 - contiguous_region: the number of contiguous regions
 - nb_internal_nodes_from_ancestor_with_updated_weight: 
 - supporting_children: The list of children levels that support the adjacency
 - predicted_edge_age_relative_to_root: 
 - predicted_edge_lca: the deepest level where this edge is identified.

Each line in the file corresponds to one ancestral / extant adjacency. Ancestral genes for which no adjacency could be identified will be listed as single column rows.

```TSV
gene1   gene2   weight  contiguous_region       nb_internal_nodes_from_ancestor_with_updated_weight     supporting_children     predicted_edge_age_relative_to_root     predicted_edge_lca
rootHOG_7046    HOG_167759      3.0     0.0     0.0     Micromonas;Ostreococcus 0.49    Mamiellales
rootHOG_7050    HOG_172696      2.0     1.0     0.0     Micromonas;Ostreococcus 0.49    Mamiellales
rootHOG_7081    HOG_169253      3.0     2.0     0.0     Micromonas;Ostreococcus 0.49    Mamiellales
rootHOG_7087    HOG_172341      2.0     3.0     0.0     Micromonas;Ostreococcus 0.49    Mamiellales
...
```

