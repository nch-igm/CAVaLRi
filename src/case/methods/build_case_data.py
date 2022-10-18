import pandas as pd
import os
import sys
import re
import json

# load local packages
sys.path.append('../..')
from config import *

def build_case_data(case):

    # Initialize data as a dict, will be casted as a json
    # Rooting the data in genes, diseases will be added to genes
    data = {
        'case_id': case.case_id,
        'genes': {}
        }

    # Obtain a set of genes to append
    gene_df = pd.read_csv(os.path.join(config['project_root'], config['gene_info']), sep = '\t').rename(columns = {'Symbol': 'gene_name'})[['GeneID', 'gene_name']]
    gene_ids = pd.Series(gene_df.GeneID, index = gene_df.gene_name).to_dict()
    variant_df = case.genotype.scored_variants.copy()
    variant_df['gene_name'] = variant_df['SNPEFF_GENE_NAME'].str[2:-2]
    variant_df['ENTREZ_ID'] = variant_df['ENTREZ_ID'].astype('int')

    # Get gene level data by grouping variant data frame
    gene_df = variant_df.copy()[['gene_name', 'ENTREZ_ID', 'MANE_REFSEQ_NUC', 'SYMBOL',
                        'total_samples', 'alt_sample_fraction']].drop_duplicates()
    for idx, row in gene_df.iterrows():
        data['genes'][row['ENTREZ_ID']] = {
            'diseases': [],
            'variants': [],
            'gene_name': row['gene_name'],
            'SYMBOL': row['SYMBOL'],
            'ENTREZ_ID': row['ENTREZ_ID'],
            'MANE_REFSEQ_NUC': row['MANE_REFSEQ_NUC'],
            'total_samples': int(row['total_samples']),
            'alt_sample_fraction': round(row['alt_sample_fraction'], 5)
            }

    # Add each variant row to the case data
    for idx, row in variant_df.iterrows():
        variant_data = {
            'CHROM': row['CHROM'],
            'POS': row['POS'],
            'REF': row['REF'],
            'ALT': row['ALT'],
            'SNPEFF_INTERPRETED_SIMPLE_LOCATION': row['SNPEFF_INTERPRETED_SIMPLE_LOCATION'],
            'SNPEFF_ANNOTATION': row['SNPEFF_ANNOTATION'],
            'SNPEFF_ANNOTATION_IMPACT': row['SNPEFF_ANNOTATION_IMPACT'],
            'SNPEFF_INTERPRETED_EFFECT': row['SNPEFF_INTERPRETED_EFFECT'],
            'SNPEFF_INTERPRETED_HGVS_FORMAT': row['SNPEFF_INTERPRETED_HGVS_FORMAT'],
            'SNPEFF_HGVS_C': row['SNPEFF_HGVS_C'],
            'SNPEFF_HGVS_P': row['SNPEFF_HGVS_P'],
            'gnomad_ex_faf95_popmax': row['gnomad_ex_faf95_popmax'],
            'SNPdogg': row['yprob'],
            'primateAI_score': row['dbnsfp_phastcons17way_primate_rankscore'],
            'max_score': row['max_score']
        }

        # Add all sample specific data
        for c in set(variant_df.columns).intersection(set(['PROBAND', 'MOTHER', 'FATHER'])):
            sample_data = json.loads(row[c].replace("'", '"')) # Set the genotype dict strings back to dictionaries
            variant_data.update({c: sample_data})

        # Add variant to the corresponding gene key
        data['genes'][row['ENTREZ_ID']]['variants'].append(variant_data)

    # Add OMIM diseases and corresponding modes of inheritance
    moi_df = pd.read_csv(os.path.join(config['project_root'], config['moi_db']))
    moi_df['OMIM'] = moi_df['omimId'].str[5:].astype('int')
    moi_df = moi_df.drop(columns = ['omimId']).dropna(axis=0).reset_index(drop=True)
    mim2gene_df = pd.read_csv(os.path.join(config['project_root'], config['mim2gene']), sep = '\t').rename(columns = {'#MIM number': 'OMIM'})
    mim2gene_df = mim2gene_df.loc[mim2gene_df['type'] == 'phenotype']
    mim2gene_df = mim2gene_df[['OMIM', 'GeneID']].drop_duplicates()
    mim2gene_df = mim2gene_df.loc[mim2gene_df['GeneID'] != '-']
    mim2gene_df['GeneID'] = mim2gene_df['GeneID'].astype('int')
    mim2gene_df = mim2gene_df.groupby('GeneID')['OMIM'].apply(list).reset_index()

    # Iterate through each OMIM ID and add to case data
    for idx, row in mim2gene_df.iterrows():
        for gene in variant_df['ENTREZ_ID'].unique():
            if data['genes'][gene]['ENTREZ_ID'] == row['GeneID']:
                for omim in row['OMIM']:
                    moi_slice = moi_df.loc[moi_df['OMIM'] == int(omim)].reset_index(drop=True)
                    moi = moi_slice['moi'][0] if not moi_slice.empty else 'UNKNOWN'
                    data['genes'][gene]['diseases'].append({
                        'omim': omim,
                        'moi': moi
                        })

    return data