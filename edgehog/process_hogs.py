#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from ete3 import Tree
import pyham
import sys


def read_tree(tree_path):
    try:
        tree = Tree(tree_path, format=1)
    except Exception as e:
        print('error: impossible to read the species tree \'%s\'' % tree_path)
        sys.exit(e)
    return tree
            

def map_hogs_onto_tree(hogs_file, tree_path, hdf5):
    if hdf5:
        use_internal_names = True
        species_resolve_mode = "OMA"
    else:
        tree = read_tree(tree_path)
        # test if root has a name to decide whether or not internal nodes must be named by pyham
        if not tree.name:
            use_internal_names = False
        else:
            try:
                float(tree.name)
                use_internal_names = False
            except:
                use_internal_names = True
        species_resolve_mode = None
    try:
        ham = pyham.Ham(tree_file = tree_path, hog_file = hogs_file, tree_format = 'newick', 
                        use_internal_name = use_internal_names, species_resolve_mode = species_resolve_mode)
    except Exception as e:
        print('error: the mapping of HOGs onto the tree failed')
        sys.exit('pyham\'s error message: %s' % e)
    return ham
    
    
# Dictionary of genes per genome extractred from the HOGs_orthoxml
def get_hogxml_entries(ham):
    genomes = ham.get_list_extant_genomes()
    entries = dict()
    protein_id_to_entry = dict()
    for genome in genomes:
        entries[genome] = list()
        protein_id_to_entry[genome] = dict()
        for gene in genome.genes:
            gene_dict = gene.get_dict_xref()
            entries[genome].append(gene_dict)
            if not 'protId' in gene_dict:
                sys.exit('No protId tag was found for gene %s in genome %s within the orthoxml file' % (gene_dict['id'], genome.name))
            protein_id_to_entry[genome][gene_dict['protId']] = gene_dict['id']
    return entries, protein_id_to_entry

