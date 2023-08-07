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
                     'frameshift_variant','stopgain','stop_lost','startloss',
                     'stoploss','frameshift_insertion','frameshift_deletion']
    snv_variants = ['missense_variant','nonsynonymous_SNV','unknown']
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
    splicing_df = var_df[(var_df['func'].isin(splice_variants)) | (var_df['func'].str.contains('splicing'))]
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
        var_df['func'].isin(splice_variants) | (var_df['func'].str.contains('splicing'))
    )]

    # Limit SNVs to only the ones that look bad
    snv_df = snv_df.merge(genotype.snv_annotations, how = 'left', on = ['CHROM','POS','REF','ALT']).fillna(0).rename(columns = {'snv_score':'score'})
    bad_snv_df = snv_df[snv_df['score'] >= config['snv_threshold']]

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
    pathogenic_df = var_df.merge(pathogenic_df, on = index_cols[:-1])
    pathogenic_max_df = pathogenic_df.groupby(index_cols[:-1])['score'].max().reset_index()
    pathogenic_df = pathogenic_df.merge(pathogenic_max_df)

    # Save variants if the top two variants in a gene average above some threshold value
    min_threshold = min([config['snv_threshold'], config['spliceai_threshold'], config['mutpredindel_threshold']])
    scored_df = pd.concat([df[index_cols] for df in [null_df, snv_df, splicing_scored_df, inframe_indel_df, clinvar_path_df] if not df.empty])
    scored_df = scored_df.drop_duplicates().reset_index(drop=True)
    scored_df = var_df.merge(scored_df, on = index_cols[:-1])
    scored_max_df = scored_df.groupby(index_cols[:-1])['score'].max().reset_index()
    scored_df = scored_df.merge(scored_max_df)
    
    # Get disease mode of inheretences
    moi_df = pd.read_csv(os.path.join(genotype.case.cohort.root_path, config['moi_db']))
    recessive_codes = set(['AR','XLR'])
    mim2gene = pd.read_csv(os.path.join(genotype.case.cohort.root_path, config['mim2gene']), sep = '\t')
    mim2gene = mim2gene.rename(columns = {'#MIM number':'omimId'})[['omimId','GeneID']].astype({'GeneID':str})
    def map_omim(row):
        return f"OMIM:{row['omimId']}"
    mim2gene['omimId'] = mim2gene.apply(map_omim, axis = 1)
    gene_df = pd.read_csv(os.path.join(genotype.case.cohort.root_path, config['gene_info']), sep = '\t')
    gene_df = gene_df[['GeneID','Symbol']].rename(columns = {'Symbol':'geneSymbol'}).astype({'GeneID':str})
    gene_disease_df = gene_df.merge(mim2gene, on = 'GeneID').astype({'GeneID':int})
    gene_disease_df = gene_disease_df.merge(moi_df, on = 'omimId')
    rd_df = gene_disease_df[gene_disease_df['moi'].str.contains('R')]
    rd_gene_ids = list(set(gene_disease_df['GeneID']))

    # Add second worst variants to the list if they can be combined with the 
    # highest scoring variant to cause a recessive disease    
    extra_variants = []
    for g in list(pathogenic_df['GENE_ID'].unique()):
        
        # Limit to only variants for gene g
        scored_g_df = scored_df[scored_df['GENE_ID'] == g]
        
        # Get second highest values if the gene is associated with a recessive condition
        sorted_scores = sorted(scored_g_df['score'].to_list(), reverse=True)
        highest_score = scored_g_df['score'].max()
        if len(scored_g_df.index) > 1 and g in rd_gene_ids:

            second_worst_df = scored_g_df.sort_values('score', ascending = False)
            second_worst_df = second_worst_df.reset_index(drop=True).iloc[:2,:]
            
            for idx, row in second_worst_df.iterrows():
                if (row['score'] + highest_score) / 2 >= min_threshold and row['score'] > 0.05:
                    extra_variants.append(row)
        
            if len(extra_variants) > 0:
                extra_variant_df = pd.concat(extra_variants, axis = 1).T.drop_duplicates().reset_index(drop=True)
                pathogenic_df = pd.concat([pathogenic_df,extra_variant_df]).drop_duplicates().reset_index(drop=True)
        

    pathogenic_dir = '/igm/home/rsrxs003/rnb/output/BL-283/clinician/pathogenic_variants'
    path_df_path = os.path.join(pathogenic_dir, f"{genotype.case.case_id}.csv")
    pathogenic_df.to_csv(path_df_path, index = False)

    all_dir = '/igm/home/rsrxs003/rnb/output/BL-283/clinician/all_variants'
    all_df_path = os.path.join(all_dir, f"{genotype.case.case_id}.csv")
    scored_df.to_csv(all_df_path, index = False)

    var_dir = '/igm/home/rsrxs003/rnb/output/BL-283/clinician/variants'
    var_df_path = os.path.join(var_dir, f"{genotype.case.case_id}.csv")
    var_df.to_csv(var_df_path, index = False)
    
    return pathogenic_df
