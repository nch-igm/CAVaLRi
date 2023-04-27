import sys
import os
import pandas as pd
import json
import math
from scipy.stats import poisson

def filter_pathogenic_variants(genotype):
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
    null_df = var_df[var_df['exon_func'].isin(null_variants)]
    null_df['score'] = 1
    snv_df = var_df[var_df['exon_func'].isin(snv_variants)]
    inframe_indel_df = var_df[var_df['exon_func'].isin(inframe_indel_variants)]
    splicing_df = var_df[var_df['func'].isin(splice_variants)]
    clinvar_vus_df = var_df[var_df['clinvar_vus_sig']]
    clinvar_path_df = var_df[var_df['clinvar_path_sig']]
    clinvar_path_df['score'] = 1
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
    snv_df = snv_df.merge(genotype.snpdogg_annotations, how = 'left', on = ['CHROM','POS','REF','ALT']).fillna(0).rename(columns = {'snpdogg_score':'score'})
    bad_snv_df = snv_df[snv_df['score'] >= config['snpdogg_threshold']]

    # Limit splicing variants to only the ones that look bad
    splicing_scored_df = splicing_df.merge(genotype.spliceai_annotations, how = 'left', on = ['CHROM','POS','REF','ALT']).fillna(0).rename(columns = {'spliceai_score':'score'})
    bad_splicing_df = splicing_scored_df[splicing_scored_df['score'] >= config['spliceai_threshold']]


    # Limit inframe indels to only those that look bad
    inframe_indel_df = inframe_indel_df.merge(genotype.mutpredindel_annotations, how = 'left', on = ['CHROM','POS','REF','ALT']).fillna(0).rename(columns = {'mutpred_score':'score'})
    bad_inframe_indel_df = inframe_indel_df[inframe_indel_df['score'] >= config['mutpredindel_threshold']]

    # Potentially pathogenic variant counts by category
    index_cols = ['CHROM','POS','REF','ALT','score']
    pathogenic_df = pd.concat([df[index_cols] for df in [null_df, bad_snv_df, bad_splicing_df, bad_inframe_indel_df, clinvar_path_df] if not df.empty])
    pathogenic_df = pathogenic_df.drop_duplicates().reset_index(drop=True)

    pathogenic_dir = '/igm/home/rsrxs003/rnb/output/BL-283/clinician/pathogenic_variants'
    path_df_path = os.path.join(pathogenic_dir, f"{genotype.case.case_id}.csv")
    var_df.merge(pathogenic_df, on = index_cols[:-1]).to_csv(path_df_path, index = False)

    scored_df = pd.concat([df[index_cols] for df in [null_df, snv_df, splicing_scored_df, inframe_indel_df, clinvar_path_df] if not df.empty])
    all_dir = '/igm/home/rsrxs003/rnb/output/BL-283/clinician/all_variants'
    all_df_path = os.path.join(all_dir, f"{genotype.case.case_id}.csv")
    var_df.merge(scored_df, on = index_cols[:-1]).to_csv(all_df_path, index = False)
    
    return var_df.merge(pathogenic_df, on = index_cols[:-1])
