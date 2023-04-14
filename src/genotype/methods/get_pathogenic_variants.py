import sys
import os
import pandas as pd
import json
import math
from scipy.stats import poisson


def get_pathogenic_variants(genotype):
    """
    input:
        genotype -- CAVaLRi Genotype type object 
    output:
        dictionary keyed by gene symbol with a float value corresponding to the 
        log_10 scaled genotype likelihood ratio
    """

    # Read in the configuration settings of the Cohort object
    config = genotype.case.cohort.config

    # Intialize result dict
    res = {}

    # Define a refgene functional list to identify null variants according to
    # https://annovar.openbioinformatics.org/en/latest/user-guide/gene/
    null_variants = ['frameshift_elongation','frameshift_truncation',
                     'frameshift_variant','stopgain','stop_lost',
                     'frameshift_insertion','frameshift_deletion']
    snv_variants = ['missense_variant','nonsynonymous_SNV']
    inframe_indel_variants = ['nonframeshift_insertion','nonframeshift_deletion',
                          'inframe_insertion','inframe_deletion',
                          'inframe_variant']
    splice_variants = ['splicing']

    # Define clinical interpretation annotation categories
    clinvar_sig = {
        'pathogenic':['Pathogenic','Likely_pathogenic','Pathogenic/Likely_pathogenic'],
        'benign_significance':['Benign','Likely_benign','Benign/Likely_benign'],
        'unknown_significance':['Conflicting_interpretations_of_pathogenicity']
    }

    # For each disease, access gene data and calculate genotype likelihood ratio

    var_df = genotype.variants.copy()
    # var_df = var_df[var_df['GENE'] == g]

    # Break variants into categories based on exon function
    def get_func(row):
        return json.loads(row['INFO'])['Func.refGene'][0]
    
    def get_exon_func(row):
        return json.loads(row['INFO'])['ExonicFunc.refGene'][0]
    
    def get_clinvar_vus_sig(row, clinvar_sig):
        info = json.loads(row['INFO'])['CLNSIG'][0]
        res = False
        if info == None:
            return res
        for cs in info.split('|'):
            if cs in clinvar_sig['unknown_significance']:
                res = True
        return res
    
    def get_clinvar_path_sig(row, clinvar_sig):
        info = json.loads(row['INFO'])['CLNSIG'][0]
        res = False
        if info == None:
            return res
        for cs in info.split('|'):
            if cs in clinvar_sig['pathogenic']:
                res = True
        return res
    
    var_df['func'] = var_df.apply(get_func, axis = 1)
    var_df['exon_func'] = var_df.apply(get_exon_func, axis = 1)
    var_df['clinvar_vus_sig'] = var_df.apply(get_clinvar_vus_sig, clinvar_sig = clinvar_sig, axis = 1)
    var_df['clinvar_path_sig'] = var_df.apply(get_clinvar_path_sig, clinvar_sig = clinvar_sig, axis = 1)
    var_df.to_csv('/Users/rsrxs003/projects/CAVaLRi_/all_var.csv', index = False)
    null_df = var_df[var_df['exon_func'].isin(null_variants)]
    snv_df = var_df[var_df['exon_func'].isin(snv_variants)]
    inframe_indel_df = var_df[var_df['exon_func'].isin(inframe_indel_variants)]
    splicing_df = var_df[var_df['func'].isin(splice_variants)]
    clinvar_vus_df = var_df[var_df['clinvar_vus_sig']]
    clinvar_path_df = var_df[var_df['clinvar_path_sig']]
    other_df = var_df[~(
        var_df['exon_func'].isin(null_variants)
            |
        var_df['exon_func'].isin(snv_variants)
            |
        var_df['exon_func'].isin(inframe_indel_variants)
            |
        var_df['func'].isin(splice_variants)
    )]

    # Limit SNVs to only the ones that look bad
    threshold = 0.8
    def get_snv_score(row):
        # return row['SNPdogg']
        return 1
    snv_df['score'] = snv_df.apply(get_snv_score, axis = 1)
    bad_snv_df = snv_df[snv_df['score'] >= config['snpdogg_threshold']]

    # Limit splicing variants to only the ones that look bad
    splicing_scored_df = splicing_df.merge(genotype.spliceai_annotations, how = 'left', on = ['CHROM','POS','REF','ALT'])
    bad_splicing_df = splicing_scored_df[splicing_scored_df['spliceai_score'] >= config['spliceai_threshold']]


    # Limit inframe indels to only those that look bad
    def get_inframe_sig(row):
        # return json.loads(row['INFO'])['inframe_sig'][0]
        return 0.2
    inframe_indel_df['inframe_score'] = inframe_indel_df.apply(get_inframe_sig, axis = 1)
    bad_inframe_indel_df = inframe_indel_df[inframe_indel_df['inframe_score'] >= config['mutpredindel_threshold']]

    # Potentially pathogenic variant counts by category
    index_cols = ['CHROM','POS','REF','ALT']
    pathogenic_df = pd.concat([
        null_df[index_cols],
        bad_snv_df[index_cols],
        bad_splicing_df[index_cols],
        bad_inframe_indel_df[index_cols],
        clinvar_path_df[index_cols]
    ]).drop_duplicates().reset_index(drop=True)
    
    return var_df.merge(pathogenic_df, on = index_cols)
