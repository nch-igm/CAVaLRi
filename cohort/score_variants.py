import sys
import os
import pandas as pd
import numpy as np

sys.path.append('../../workflow')
from config import *

def get_score(row):
    
    # Intialize result
    res = []

    # Iterate through each score
    for score in row['score']:
        if score == 'SNPDogg':
            res.append(max([0.2, row['yprob']]))
        elif score == 'SpliceAI':
            res.append(0.5)
        elif score == 'PrimateAI, SURF':
            res.append(0.2)
        elif score == 'MutPredAI': 
            res.append(0.4)
        else:
            res.append(float(score))
        
    return max(res)


def score_variants(cohort):

    # Expand variants with multiple snpeff annotations
    s = cohort.annotated_variants['SNPEFF_ANNOTATION'].str.split('&').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'SNPEFF_ANNOTATION'
    res_df =  cohort.annotated_variants.drop(columns = ['SNPEFF_ANNOTATION']).join(s)

    # Append scores
    var_scores = pd.read_csv(os.path.join(config['project_root'], config['genotype_scores']))
    res_df = res_df.merge(var_scores[['seq_ontology', 'score']], left_on = 'SNPEFF_ANNOTATION', right_on = 'seq_ontology', how = 'left')

    # Change dictionary columns to string columns
    for c in set(res_df.columns).intersection(set(['PROBAND', 'MOTHER', 'FATHER'])):
        res_df[c] = res_df[c].astype('str')

    # Group scores
    res_df = res_df.fillna(False).groupby([item for item in res_df.columns if item not in ['score']])['score'].apply(list).reset_index(name = 'score')
    res_df['max_score'] = res_df.apply(get_score, axis = 1)

    res_df.reset_index(drop=True).to_csv('hi_df.csv', index=False)

    return res_df.drop(columns = ['seq_ontology', 'score']).drop_duplicates().reset_index(drop=True)
