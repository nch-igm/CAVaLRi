import os
import pandas as pd
import json
import pickle

def build_case_data(case):

    # Intialize case data
    case_data = {'genes': {k:{} for k in case.genotype.variants['GENE'].unique()}}
    
    # Read in gene_disease dataframe
    gene_path = '/Users/rsrxs003/projects/CAVaLRi_/data/Homo_sapiens_gene_info.tsv'
    gene_df = pd.read_csv(gene_path, sep = '\t').astype({'GeneID': str})
    gene_disease_path = '/Users/rsrxs003/projects/LIRICAL/data/mim2gene_medgen'
    gene_disease_df = pd.read_csv(gene_disease_path, sep = '\t').rename(columns = {'#MIM number':'OMIM'})
    # with open('/Users/rsrxs003/projects/CAVaLRi_/example/check.txt','w') as f:
        # print(gene_disease_df.dtypes, file = f)
        # print(gene_df.dtypes, file = f)
        # case_data
    gene_disease_df = gene_disease_df.merge(gene_df[['GeneID','Symbol']], on = 'GeneID')
    gene_disease_df = gene_disease_df[gene_disease_df['type'] == 'phenotype']

    for k,v in case_data['genes'].items():

        omim_ids = gene_disease_df[gene_disease_df['Symbol'] == k]['OMIM'].to_list()
        for oi in omim_ids:
            v[oi] = {}
    
    case_data = {'genes': {k:v for k,v in case_data['genes'].items() if len(case_data['genes'][k]) > 0}}
    with open('/Users/rsrxs003/projects/CAVaLRi_/example/case_data.json', 'w') as f:
        json.dump(case_data, f, indent = 4)
    
    case.case_data_ = case_data
    with open('/Users/rsrxs003/projects/CAVaLRi_/example/case_earlier_.pickle', 'wb') as f:
        pickle.dump(case, f)

    return case_data
    