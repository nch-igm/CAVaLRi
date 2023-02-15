import pandas as pd
import os
import sys
import re
import json


def read_phenotypes(case):

    # Intialize phenotype dictionary
    res = {}

    # Read in HPO IC
    ic_df = pd.read_csv(os.path.join(case.cohort.root_path,'data','HPO_with_gene_IC.csv'))

    # Read in phenotypes stored as csvs
    pheno_df = pd.read_csv(case.phenotype.phenotype_path).merge(ic_df)
    
    pheno_df = pheno_df.sort_values('IC', ascending = False).reset_index(drop = True).reset_index().rename(columns = {'index': 'rank'})
    pheno_df['rank'] = pheno_df['rank'] + 1

    for idx, row in pheno_df.iterrows():

        # Add HPO ID to the result
        res[row['HPO ID']] = row['rank']
        
    return res, pheno_df
