#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import itertools
import numpy
import pandas
import sys


# split the column attributes of the gff using the ";" field separator
# string before "=" becomes the column name and string after "=" becomes the variable
# returns the attributes as a pandas dataframe
def gff_attributes_to_columns(df):
    attributes_df = pandas.DataFrame(df['attributes'].apply(
        lambda attributes: dict([key_value_pair.split(sep='=', maxsplit=1) for key_value_pair in attributes.strip(';').split(';')])))
    attributes_df.columns = ['at_dic']
    attributes_df['at_dic_keys'] = attributes_df['at_dic'].apply(lambda at_dic: list(at_dic.keys()))
    merged_attributes_list = list(itertools.chain.from_iterable(attributes_df['at_dic_keys']))
    nonredundant_list = sorted(list(set(merged_attributes_list)))
    for atr in nonredundant_list:
        df[atr] = attributes_df["at_dic"].apply(lambda at_dic: at_dic.get(atr))
    return df


# Given a sorted dataframe of annotations limited to (pseudo)genes and cds,
# this function associates the parent (pseudo)gene to each CDS
# this is done based on the indices of (pseudo)genes in the df:
# if one or more CDS(s) appear(s) just after a (pseudo)gene line, 
# then the cds is/are considered part of the gene until another gene occurs
#
# returns a filtered_df (with non_coding genes and pseudo_genes removed)
# with the novel columns 'gene_index', 'parent_feature', 'position'
def get_gene_position_and_cds_children(df):
    df['feature_index'] = df.index
    features = set(df['type'])
    if 'GENE' in features or 'Gene' in features or 'gene' in features:
        only_cds = False
        gene_df = df.loc[df['type'].isin(['GENE', 'Gene', 'gene', 'PSEUDOGENE', 'Pseudogene', 'pseudogene'])]
        # gene_id is a vector of the same length as nb of rows in df.
        # a specific 'gene_id' is repeated until another gene or pseudogene is encountered in the df
        # when compared side by side to protein_id, this allows to tell which gene is the parent of each cds
        gene_id = numpy.zeros((df.shape[0], ), dtype=int)
        idx_old_gene = 0
        for idx_new_gene in gene_df['feature_index']:
            gene_id[(idx_old_gene):idx_new_gene+1] = idx_old_gene
            idx_old_gene = idx_new_gene
        gene_id[(idx_old_gene):] = idx_old_gene
        df['gene_index'] = gene_id
        df['parent_feature'] = df['type'][gene_id].tolist()
        # retain only coding genes (those encompassing CDSs)
        df = df.loc[df['gene_index'].isin(df.loc[df['type'].isin(['CDS', 'Cds', 'cds'])]['gene_index'])]
        # remove pseudogenes
        df = df.loc[df['parent_feature'].isin(['GENE', 'Gene', 'gene'])].reset_index()
        # reassign gene index, now that non-coding genes and pseudogenes have been filtered out 
        uniq_genes = set(df['gene_index'])
        g2i = dict(zip(uniq_genes, range(0, len(uniq_genes))))
        df['gene_index'] = [g2i[g] for g in df['gene_index']]
    else:
        only_cds = True # this is the case for instance in the output of prodigal
        df['gene_index'] = df['feature_index']
    # get position of each gene on each contig (genomic_accession) present in the gff
    # (in fact the code below is not necessary since only the order on a contig matter.
    # whether first gene has an index of 0 or 1000 is essentially the same)
    # contig_first_gene = df.drop_duplicates(subset='genomic_accession')['gene_index']
    # id_old_contig_first_gene = 0
    # row_old_contig_first_gene = 0
    # position_subtrahend = numpy.zeros((df.shape[0], ), dtype=int)
    # for i in range(1, len(contig_first_gene)):
    #     id_new_contig_first_gene = contig_first_gene.iloc[i]
    #     row_new_contig_first_gene = contig_first_gene.index[i]
    #     position_subtrahend[row_old_contig_first_gene:row_new_contig_first_gene] = id_old_contig_first_gene
    #     id_old_contig_first_gene = id_new_contig_first_gene
    #     row_old_contig_first_gene = row_new_contig_first_gene
    # position_subtrahend[row_old_contig_first_gene:] = id_old_contig_first_gene
    # df['position'] = df['gene_index'] - position_subtrahend
    df['position'] = df['gene_index']
    return df, only_cds
     

# make sure the main_key is 'protein_id' for proteins and 'gene_id' for genes
def make_key(df, main_key, alternative_keys):
    if not main_key in df:
        for key in alternative_keys:
            if key in df:
                df[main_key] = df[key]
                break
    return df
            

# turns a gff into a pandas dataframe,
# modify the df to include key informations such as the parent gene of each cds 
# and the position of each gene on each contig
# Eventually, the data are integrated into python dictionaries ('gene_dict' and 'protein_dict')
def gff_to_dict(gff):
    gff_column_names = ['genomic_accession', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    useful_column_names = ['genomic_accession', 'type', 'start', 'end', 'strand', 'attributes']
    column_types = {'genomic_accession': str, 'type': str, 'start': int, 'end': int, 'strand': str, 'attributes': str}
    try:
        df = pandas.read_csv(gff, sep='\t', comment='#', header=None, names=gff_column_names, usecols=useful_column_names, dtype=column_types)
    except Exception as e:
        print('EXIT: \'%s\' does not comply with gff format requirements' % gff)
        print('Reading the gff exited with the following error:')
        sys.exit(e)
    # retain only CDS and genes
    df = df.loc[df['type'].isin(['CDS', 'Cds', 'cds', 'GENE', 'Gene', 'gene', 'PSEUDOGENE', 'Pseudogene', 'pseudogene'])]
    # check there are cds or genes in the gff
    if df.empty:
        sys.exit('error: column \'type\' (fields n°3) of the input gff never corresponds to the \'CDS\' or \'gene\' string')
    
    # turn attributes string into columns
    # df['protein_id'] = df['attributes'].apply(
    #     lambda attributes: (attributes + ";protein_id=").strip(';').split(sep = "protein_id=", maxsplit=1)[1].split(";")[0])
    df = gff_attributes_to_columns(df)
    
    df.reset_index(inplace=True)
    df, only_cds = get_gene_position_and_cds_children(df)
    df_proteins = df.loc[df['type'].isin(['CDS', 'Cds', 'cds'])]
    if not 'protein_id' in df_proteins:
        if 'ID' in df_proteins:
            temp_id = df_proteins['ID'] 
        elif 'id' in df:
            temp_id = df_proteins['id']
        else:
            sys.exit('Neither the tag \'protein_id\' nor \'ID\' are present in the column \'attributes\' of the input gff: \'%s\'' % gff)
        df_proteins['protein_id'] = df['genomic_accession'] + temp_id.apply(lambda x: '_' + x.split('_')[1]) 
    if only_cds:
        df['gene_id'] = df['protein_id']
        df_genes = df
    else:
        df = make_key(df, 'gene_id', ['GeneID', 'geneID', 'gene_id', 'ID', 'Id', 'id', 'Name', 'name', 'gene_index'])
        df_genes = df.loc[df['type'].isin(['GENE', 'Gene', 'gene'])]   
    df_proteins = df_proteins.dropna(subset=['protein_id']).drop_duplicates('protein_id', keep='first')
    gene_dict = df_genes[['genomic_accession', 'gene_index', 'gene_id', 'position', 'start', 'end', 'strand']].set_index(['genomic_accession', 'gene_index']).T.to_dict()
    protein_dict = df_proteins[['genomic_accession', 'protein_id', 'start', 'end', 'strand', 'gene_index', 'position']].set_index('protein_id').T.to_dict()
    
    # position_to_gene = df_genes[['genomic_accession', 'position', 'gene_index']].set_index(['genomic_accession', 'position']).T.to_dict()
    return gene_dict, protein_dict 
 
    