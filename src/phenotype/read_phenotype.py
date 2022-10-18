import pandas as pd
import os
import sys
import re
import json

# load local packages
sys.path.append('../..')
from config import *


def read_phenotypes(case):

    # Intialize phenotype dictionary
    res = {}

    # Read in phenotypes stored as csvs
    pheno_df = pd.read_csv(case.phenotype.phenotype_path).reset_index().rename(columns = {'index': 'rank'})
    pheno_df['rank'] = pheno_df['rank'] + 1

    for idx, row in pheno_df.iterrows():
        
        # Only keep the top n phenotypes (where n = phenotype set length upper limit)
        if row['rank'] > config['hpo_total_upper_bound']:
            break

        # Add HPO ID to the result
        res[row['HPO ID']] = row['rank']
        
    
    return res


# def read_phenotype(case):
#     phenotypes = 'hi'
#     omim_path = os.path.join(config['project_root'], config['hpoa'])
#     with open (omim_path, 'r') as f:
#         for line in f:
#             if re.search('^#DatabaseID', line):
#                 columns = line[1:].replace('\n', '').split('\t')

#     omim_df = pd.read_csv(omim_path, comment = '#', sep = '\t')
#     omim_df.columns = columns