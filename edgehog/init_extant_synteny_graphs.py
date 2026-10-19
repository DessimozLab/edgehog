#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import networkx as nx
# from multiprocessing import Pool
from edgehog.read_gff import gff_to_dict


# for each genome, identify <gff_directory>/<genome_name>*.gff(3), optionally compressed (.gz, .bz2, .xz, .zip),
# and exit if there is none or several; other files such as checksums or indexes (.md5, .tbi) are ignored
# a file is attributed to the longest genome name it starts with (HUMAN2.gff belongs to HUMAN2, not HUMAN),
# and <genome_name>.gff(3) is preferred over other files starting with the genome name
# returns a genome_to_gff dictionary
def get_gffs_of_extant_genomes(gff_directory, genomes):
    files = sorted(os.listdir(gff_directory))
    patterns = {g: re.compile(re.escape(g.name) + r'(.*?)\.gff3?(\.(gz|bz2|xz|zip))?$') for g in genomes}
    file_owner = dict()
    for f in files:
        owners = [g for g in genomes if patterns[g].match(f)]
        if owners:
            file_owner[f] = max(owners, key=lambda g: len(g.name))
    genome_to_gff_dict = dict()
    for g in genomes:
        candidates = [f for f in files if file_owner.get(f) is g]
        exact = [f for f in candidates if patterns[g].match(f).group(1) == '']
        if len(exact) == 1:
            match = exact[0]
        elif len(exact) == 0 and len(candidates) == 1:
            match = candidates[0]
        elif len(candidates) == 0:
            sys.exit('\'%s.gff\' does not exist in the \'%s\' gff_directory' % (g.name, gff_directory))
        else:
            sys.exit('several gffs found for genome \'%s\' in the \'%s\' gff_directory: %s'
                     % (g.name, gff_directory, ', '.join(exact or candidates)))
        genome_to_gff_dict[g] = os.path.join(os.path.abspath(gff_directory), match)
    return genome_to_gff_dict


# the protein_id in hogs_xml is basically the header of the protein sequences in the fasta given as input to OMA standalone
# this protein_id is then likely to be longer string than the true protein_id written in the gff
# this function modifies the protein_id in the dictionary derived from the hogs_xml 
# in order to properly intersect it with the dictionaries derived from the gff
def intersect_hogxml_and_gff(genome, hogxml_entries, protein_id_to_hogxml_entry, protein_dict, gene_dict):
    hogxml_protein_ids = list(protein_id_to_hogxml_entry.keys())
    if hogxml_protein_ids[0] in protein_dict:
        real_protein_ids = hogxml_protein_ids
    elif hogxml_protein_ids[0].split(' ')[0] in protein_dict:
        real_protein_ids = [p.split(' ')[0] for p in hogxml_protein_ids]
    elif hogxml_protein_ids[0].split(',')[0] in protein_dict:
        real_protein_ids = [p.split(',')[0] for p in hogxml_protein_ids]
    else:
        sys.exit('error: unable to properly intersect hogs_xml with gff for genome %s' % genome.name)
    hogxml_protein_id_to_real_protein_id = dict(zip(hogxml_protein_ids, real_protein_ids))
    w = 0
    for p in hogxml_protein_ids:
        if w < 50:
            entry = protein_id_to_hogxml_entry[p]
            real_protein_id = hogxml_protein_id_to_real_protein_id[p]
            if real_protein_id in protein_dict:
                protein_attributes = protein_dict[real_protein_id]
                genomic_accession = protein_attributes['genomic_accession']
                gene_index = protein_attributes['gene_index']
                if not 'hogxml_entries' in gene_dict[(genomic_accession, gene_index)]:
                    gene_dict[(genomic_accession, gene_index)]['hogxml_entries'] = list()
                gene_dict[(genomic_accession, gene_index)]['hogxml_entries'].append(entry)
            else:
                print('unable to find %s in the gff of genome %s' % (p, genome.name))
                w += 1;
        else:
            sys.exit('error: unable to properly intersect hogs_xml with gff for genome %s' % genome.name)
    return gene_dict


# if a gene encodes several proteins (isoforms), 
# choose the most represented HOG among these proteins as the reference for the gene
def get_representative_hog_for_gene(ham, hogxml_entries):
    occurrence = dict()
    best_hog_occurrence = 0
    best_hog = None
    for entry in hogxml_entries:
        gene = ham.get_gene_by_id(entry)
        if gene.parent is not None: 
            hog = ham.get_hog_by_gene(gene)
            if hog:
                if hog in occurrence:
                    occurrence[hog] += 1
                else:
                    occurrence[hog] = 1
                if occurrence[hog] > best_hog_occurrence:
                    best_hog_occurrence = occurrence[hog]
                    best_hog = gene
    if best_hog:
        return best_hog
    return gene
            

# Create a feature "synteny" that stores a synteny graph for the current extant genome
def assign_extant_synteny(genome, gene_dict, ham, orient_edges): 
    graph = nx.Graph(genome=genome)
    contiguity_dict = dict()
    # gffs are not necessarily sorted: order the genes by contig (in order of first appearance in the gff)
    # and by genomic coordinates on the contig
    contig_rank = dict()
    for contig, _ in gene_dict:
        contig_rank.setdefault(contig, len(contig_rank))
    gene_keys = sorted(gene_dict, key=lambda key: (contig_rank[key[0]], gene_dict[key]['start'], gene_dict[key]['end']))

    # find first gene with hogxml_entries, i.e. that is reported in orthoxml
    for first, gene_key in enumerate(gene_keys):
        if 'hogxml_entries' in gene_dict[gene_key]:
            break
    else:
        raise ValueError('No gene with hogxml_entries found in genome %s' % genome.name)

    old_key = gene_keys[first]
    old_contig = old_key[0]
    hogxml_entries = gene_dict[gene_keys[first]]['hogxml_entries']
    old_gene = get_representative_hog_for_gene(ham, hogxml_entries)
    graph.add_node(old_gene, **gene_dict[gene_keys[first]], contig=old_contig)
    for i in range(first+1, len(gene_keys)):
        contig = gene_keys[i][0]
        if 'hogxml_entries' in gene_dict[gene_keys[i]]:
            hogs_entries = gene_dict[gene_keys[i]]['hogxml_entries']
            gene = get_representative_hog_for_gene(ham, hogs_entries)
            graph.add_node(gene, **gene_dict[gene_keys[i]], contig=contig)
            if contig == old_contig:
                # connect genes only if they are on the same contig
                graph.add_edge(old_gene, gene, weight=1, unidirectional=0, convergent=0, divergent=0, children=[None], extant_descendants=[None])
                if orient_edges:
                    # strand of the previous gene in the orthoxml, not of a skipped gene in between
                    old_gene_strand = gene_dict[old_key]['strand']
                    gene_strand = gene_dict[gene_keys[i]]['strand']
                    try:
                        if old_gene_strand == gene_strand:
                            graph[old_gene][gene]['unidirectional'] = 1
                        elif old_gene_strand == "+":
                            graph[old_gene][gene]['convergent'] = 1
                        elif old_gene_strand == "-":
                            graph[old_gene][gene]['divergent'] = 1
                    except:
                        pass
            old_gene, old_contig, old_key = gene, contig, gene_keys[i]
    # for cc in nx.connected_components(graph):
    #     if len(cc > 1):
    #         for i in range(0, len(cc)):
    #             for j in range(i+1, len(cc)):
    #                 node_a, node_b = cc[i], cc[j]
    #                 contiguity_dict[(node_a, node_b)] = contiguity_dict[(node_b, node_a)] = 1   
    return graph, contiguity_dict


def gff_to_graph(ham, genome, gff, hogxml_entries, protein_id_to_hogxml_entry, orient_edges):
    gene_dict, protein_dict = gff_to_dict(gff)
    gene_dict = intersect_hogxml_and_gff(genome, hogxml_entries[genome], protein_id_to_hogxml_entry[genome], protein_dict, gene_dict)
    graph, contiguity_dict = assign_extant_synteny(genome, gene_dict, ham, orient_edges)
    return graph, contiguity_dict
    
    
def init_extant_graphs(cpu, ham, hogs_file, gff_directory, hogxml_entries, protein_id_to_hogxml_entry, orient_edges):
    print('###################################')
    print('Initializing synteny graphs of extant genomes ...')
    ext_genomes = ham.get_list_extant_genomes()    
    genome_to_gff_dict = get_gffs_of_extant_genomes(gff_directory, ext_genomes)
    if cpu == 1:
        graphs = list()
        contiguity_dicts = list()
        i = 1
        for genome in ext_genomes:
            print('processing extant genome %d/%d: %s' % (i, len(ext_genomes), genome.name))
            graph, contiguity_dict = gff_to_graph(ham, genome, genome_to_gff_dict[genome], hogxml_entries, protein_id_to_hogxml_entry, orient_edges)
            graphs.append(graph)
            contiguity_dicts.append(contiguity_dict)
            i += 1
    # else:
    #     n = len(ext_genomes)
    #     gffs_pool = Pool(cpu)
    #     gffs = [genome_to_gff_dict[genome] for genome in ext_genomes]
    #     input_variables = list(zip([ham]*n, ext_genomes, gffs, [hogxml_entries]*n, [protein_id_to_hogxml_entry]*n))
    #     graphs = gffs_pool.starmap(gff_to_graph, input_variables)
        # problem correspondance graph and genome in case of parallelization
    for i in range(0, len(graphs)):
        g = graphs[i]
        contiguity_dict = contiguity_dicts[i]
        genome = g.graph['genome']
        genome.taxon.add_prop('bottom_up_synteny', g)
        genome.taxon.add_prop('contiguity_dict', contiguity_dict)  
    return ham


def init_extant_graphs_from_hdf5(ham, hdf5_file, orient_edges):
    from .hdf5 import init_extant_graphs_from_hdf5 as extant_graphs_from_hdf5
    return extant_graphs_from_hdf5(ham, hdf5_file, orient_edges)


