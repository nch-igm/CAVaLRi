import sys
import os
import pandas as pd
import yaml
import json
import argparse

sys.path.append('..')
sys.path.append('.')

from config import *

def build_case_summary(case_data):

    # Initialize result data frame
    res_df = pd.DataFrame(columns=['geneRank', 'geneNcbiId', 'geneSymbol', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR', 'hpoCount'])
    # res_df = pd.DataFrame(columns=['geneRank', 'geneNcbiId', 'geneSymbol', 'postTestProbability', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR', 'hpoCount'])

    for gene in case_data['genes']:
        for d in case_data['genes'][gene]['diseases']:
            res_df = res_df.append({
                'geneRank': d['geneRank'],
                'geneNcbiId': gene,
                'geneSymbol': case_data['genes'][gene]['SYMBOL'],
                # 'postTestProbability': d['postTestProbability'],
                'phenoLR': d['phenoLR_log10'],
                'geneLR': case_data['genes'][gene]['genotypeLR_log10'],
                'moiLR': d['moiLR_log10'],
                'compositeLR': d['compositeLR_log10'],
                'hpoCount': d['phenoCount']
            }, ignore_index = True)

    # Scale likelihood ratios
    res_df['phenoLR'] = res_df['phenoLR'] * config['phenoLR_scalar']
    res_df['geneLR'] = res_df['geneLR'] * config['geneLR_scalar']
    res_df['moiLR'] = res_df['moiLR'] * config['moiLR_scalar']
    res_df['compositeLR'] = res_df['phenoLR'] + res_df['geneLR'] + res_df['moiLR']

    # Get highest compositeLR for each disease
    compositLR_max_df = res_df[['geneSymbol', 'compositeLR']].groupby('geneSymbol').max().reset_index().rename(columns=({'compositeLR': 'max_compositeLR'}))
    res_df = res_df.merge(compositLR_max_df)
    res_df = res_df.groupby(['geneRank', 'geneNcbiId', 'geneSymbol', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR', 'max_compositeLR']).first().reset_index()
    res_df = res_df.loc[res_df['compositeLR'] == res_df['max_compositeLR']].drop(columns=['max_compositeLR']).drop_duplicates().sort_values('geneRank').reset_index(drop=True)

    return res_df

