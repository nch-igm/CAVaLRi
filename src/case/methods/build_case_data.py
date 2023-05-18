import os
import pandas as pd
import json
import pickle


def build_case_data(case):

    config = case.cohort.config
    root_path = case.cohort.root_path

    # Intialize case data
    case_data = {'genes': {k:{} for k in case.genotype.pathogenic_variants['GENE'].unique()}}

    with open('/igm/home/rsrxs003/rnb/notebooks/BL-293/cd.json', 'w') as f:
        json.dump(case_data, f, indent = 4)
    
    # Read in gene_disease dataframe
    gene_path = os.path.join(root_path, config['gene_info'])
    gene_disease_path = os.path.join(root_path, config['mim2gene'])

    gene_df = pd.read_csv(gene_path, sep = '\t').astype({'GeneID': str})
    gene_disease_df = pd.read_csv(gene_disease_path, sep = '\t').rename(columns = {'#MIM number':'OMIM'})
    gene_disease_df = gene_disease_df.merge(gene_df[['GeneID','Symbol']], on = 'GeneID')
    gene_disease_df = gene_disease_df[gene_disease_df['type'] == 'phenotype']

    # Remove disease without an annotated mode of inheretence
    moi_df = pd.read_csv(os.path.join(root_path, config['moi_db']))
    moi_diseases = [int(moi[moi.find(':')+1:]) for moi in list(set(moi_df['omimId']))]
    gene_disease_df = gene_disease_df[gene_disease_df['OMIM'].isin(moi_diseases)]
    gene_disease_df.to_csv('/igm/home/rsrxs003/rnb/notebooks/BL-293/cd.csv', index=False)

    for k,v in case_data['genes'].items():

        omim_ids = gene_disease_df[gene_disease_df['Symbol'] == k]['OMIM'].to_list()
        for oi in omim_ids:
            v[oi] = {}

    with open('/igm/home/rsrxs003/rnb/notebooks/BL-293/cd1.json', 'w') as f:
        json.dump(case_data, f, indent = 4)
    
    case_data = {'genes': {k:v for k,v in case_data['genes'].items() if len(case_data['genes'][k]) > 0}}

    with open('/igm/home/rsrxs003/rnb/notebooks/BL-293/cd2.json', 'w') as f:
        json.dump(case_data, f, indent = 4)

    return case_data
    