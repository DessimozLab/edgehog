# Example dataset to run edgeHOG

This folder contains an example dataset to test edgeHOG. It has been assembled by extracting for 11 Similiformes genomes from the OMA Jul2023 genomes the genes located on chromosome 12.
We therefor do not expect that all genes among all genomes have homologs in this dataset. You would usually use the genes on *all* chromosomes.

The dataset consists of the following species:

```txt
           /-CALJA
          |
          |                               /-MACFA
          |                         /Macaca
          |                        |      \-MACMU
-Simiiformes         /Cercopithecinae
          |         |              |--PAPAN
          |         |              |
          |         |               \-CHLSB
          |         |
           \Catarrhini                  /-PONAB
                    |                  |
                    |          /Hominidae        /-GORGO
                    |         |        |        |
                    |         |        |        |    /-PANPA
                    |         |         \Homininae-Pan
                     \Hominoidea                |    \-PANTR
                              |                 |
                              |                  \-HUMAN
                              |
                               \-NOMLE
```

The chromosome 12 of these genomes contains the following number of genes:

| Genome | Scientific name | Source | Nr Genes |
| -------|-----------------|--------|----------|
| CALJA | Callithrix jacchus | Ensembl 106 | 1110 |
| MACFA | Macaca fascicularis | Ensembl 94 | 722 |
| MACMU | Macaca mulatta | Ensembl 106 | 726 |
| PAPAN | Papio anubis | Ensembl 106 | 1275 |
| CHLSB | Chlorocebus sabaeus | Ensembl 77 | 706 |
| PONAB | Pongo abelii | Ensembl 106 | 1063 |
| GORGO | Gorilla gorilla gorilla | Ensembl 106 | 1121 |
| PANPA | Pan paniscus | Ensembl Main 91 | 1124 |
| PANTR | Pan troglodytes | Ensembl Main 91 | 1149 |
| HUMAN | Homo sapiens | Ensembl 102 | 1033 |
| NOMLE | Nomascus leucogenys | Ensembl 106 | 1114 |
------------------------------------------------------

We inferred HOGs using FastOMA 0.3.3 using the protein sequences available in `sequences.tgz` and the `species_tree.nwk`. 
The inferred HOGs are available in FastOMA_HOGs.orthoxml.

EdgeHOG can be run using the FastOMA_HOGs.orthoxml file together with the 
gff3 files (available in `gff3/*gff`). Again, these contain only the genes located on chromosome 12.

The following command will compute the ancestral edges with `edgeHOG`:

```bash
edgehog  --hog test_data/FastOMA_HOGs.orthoxml \
         --species_tree test_data/species_tree.nwk \
         --gff_directory test_data/gff3/ \
         --date_edges \
         --output_directory test-results
```
