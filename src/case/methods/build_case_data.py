import os
import pandas as pd
import json
import pickle


def build_case_data(case):

    config = case.cohort.config
    root_path = case.cohort.root_path

    # Get variants
    var_df = case.genotype.pathogenic_variants.copy()

    # Intialize case data
    case_data = {'genes': {k:{} for k in var_df['GENE_ID'].unique()}}
    
    # Read in gene_disease dataframe
    gene_path = os.path.join(root_path, config['gene_info'])
    gene_disease_path = os.path.join(root_path, config['mim2gene'])

    gene_df = pd.read_csv(gene_path, sep = '\t').astype({'GeneID': str})[['GeneID','Symbol']]
    gene_lookup = {row['GeneID']:row['Symbol'] for idx, row in gene_df.iterrows()}
    gene_disease_df = pd.read_csv(gene_disease_path, sep = '\t').rename(columns = {'#MIM number':'OMIM'})
    gene_disease_df = gene_disease_df.merge(gene_df, on = 'GeneID')
    gene_disease_df = gene_disease_df[gene_disease_df['type'] == 'phenotype']

    # Remove disease without an annotated mode of inheretence
    moi_df = pd.read_csv(os.path.join(root_path, config['moi_db']))
    moi_diseases = [int(moi[moi.find(':')+1:]) for moi in list(set(moi_df['omimId']))]
    gene_disease_df = gene_disease_df[gene_disease_df['OMIM'].isin(moi_diseases)]

    for k,v in case_data['genes'].items():

        omim_ids = gene_disease_df[gene_disease_df['GeneID'] == k]['OMIM'].to_list()
        for oi in omim_ids:
            v[str(oi)] = {}
    
        # Add genotype data
        g_var_df = var_df[var_df['GENE_ID'] == k].reset_index(drop=True)
        v['gene_data'] = {
            'symbol': gene_lookup[k],
            'variants': [
                {
                    'hg38_position': f"{row['CHROM']}:{row['POS']}{row['REF']}>{row['ALT']}",
                    'score': row['score'],
                    'clinvar_pathogenic': row['clinvar_path_sig'],
                    'function': row['func'],
                    'exon_function': row['exon_func']
                }
            for idx,row in g_var_df.iterrows()]
        }
    
    case_data = {'genes': {k:v for k,v in case_data['genes'].items() if len(v.keys()) > 1}}

    return case_data
    