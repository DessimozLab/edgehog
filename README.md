# edgeHOG

```edgeHOG``` is a tool to infer the gene order of each ancestor in a species phylogeny. As such, ```edgeHOG``` enables both to explore ancestral microsyntenies (local scale) and to reconstruct ancestral chromosomes (global scale). 

```edgeHOG``` relies on objects called HOGs (Hierarchical Groups of Orthologs) to model gene lineages and ancestral gene content. Basically, genes that belong to the same HOG across extant genomes are inferred to have descended from the same common ancestral gene in the common ancestor of these genomes. Accordingly, adjacencies between extant genes can be converted to edges between HOGs, which enables parsimonious ancestral gene order inferences.  

## Installation

### From PyPi using pip
`edgeHOG`  can be installed directly from pypi using pip. The command is the following:

```bash
pip install edgehog
```

### From sources
```edgeHOG``` was built and tested with python 3.6. To set up ```edgeHOG``` on your local machine, please follow the instructions below

```bash
git clone https://github.com/dessimozlab/edgehog.git
cd edgehog
poetry install
```

## Usage

> [!CAUTION]
> This branch `oma` contains adaptions to integrate the results directly into a OmaServer.h5 file.
> You most likely should switch to the main branch.

```
usage: edgehog [-h] [--version] [--output_directory OUTPUT_DIRECTORY]
               --species_tree SPECIES_TREE --hogs HOGS
               [--gff_directory GFF_DIRECTORY] [--hdf5 HDF5] [--date_edges]
               [--max_gaps MAX_GAPS]

edgehog is a software tool that infers an ancestral synteny graph at each
internal node of an input species phylogenetic tree

optional arguments:
  -h, --help            show this help message and exit
  --version             print version number and exit
  --output_directory OUTPUT_DIRECTORY
                        path to output directory (default is current
                        directory)
  --species_tree SPECIES_TREE
                        path to species/genomes phylogenetic tree (newick
                        format)
  --hogs HOGS           path to the HierarchicalGroups.orthoxml file in which
                        HOGs are stored
  --gff_directory GFF_DIRECTORY
                        path to directory with the gffs of extant genomes
                        (each gff file must be named according to the name of
                        an extant genome / leaf on the species tree)
  --hdf5 HDF5           path to the hdf5 file (alternative to gff_directory to
                        run edgeHOG on the OMA database)
  --date_edges          whether the age of edges in extant species should be
                        predicted

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

