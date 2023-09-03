import pandas as pd
import os
import sys
import re
import json


def read_phenotypes(case):

    config = case.cohort.config
    root_path = case.cohort.root_path

    # Intialize phenotype dictionary
    res = {}

    # Read in default HPO scores
    ic_df = pd.read_csv(os.path.join(root_path, config['pheno_score_source']))

    # Read in phenotypes stored as csvs with custom IC score, if available
    pheno_df = pd.read_csv(case.phenotype.phenotype_path)
    if len(pheno_df.columns) > 1:

        # Use the provided IC
        pheno_df = pheno_df.iloc[:,:2]
        pheno_df.columns = ['HPO ID','IC']
        pheno_df = pheno_df.merge(ic_df[['HPO ID','term_name']])

    else:
        
        # Read in HPO IC (Phrank score)
        pheno_df.columns = ['HPO ID']
        pheno_df = pheno_df.merge(ic_df)
        

    pheno_df = pheno_df.sort_values('IC', ascending = False).reset_index(drop = True).reset_index().rename(columns = {'index': 'rank'})
    pheno_df['rank'] = pheno_df['rank'] + 1

    for idx, row in pheno_df.head(config['hpo_upper_bound']).iterrows():
    # for idx, row in pheno_df.iterrows():

        # Add HPO ID to the result
        res[row['HPO ID']] = {
            'IC': row['IC'],
            'name': row['term_name'],
            'rank': row['rank']
        }
        
    return res, pheno_df
