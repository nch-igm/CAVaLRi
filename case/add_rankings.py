import os
import sys
import pandas as pd
import json
import argparse

# Add package locations to sys.path
sys.path.append('..')
sys.path.append('.')
from config import *


def add_rankings(case_data):

    # Intialize rank data frame
    rank_df = pd.DataFrame(columns=['omimId', 'geneId', 'compositeLR_log10'])

    for gene in case_data['genes']:
         for d in case_data['genes'][gene]['diseases']:
            rank_df = rank_df.append({
                'omimId': d['omim'],
                'geneId': gene,
                'compositeLR_log10': d['compositeLR_log10']
            }, ignore_index = True)

    # Build disease and gene ranks
    gene_rank_df = rank_df[['geneId', 'compositeLR_log10']].groupby('geneId').max().reset_index().sort_values(by='compositeLR_log10', ascending=False)
    gene_rank_df = gene_rank_df.reset_index(drop=True).reset_index().rename(columns=({'index': 'geneRank'})).drop(columns=['compositeLR_log10'])
    disease_rank_df = rank_df[['omimId', 'compositeLR_log10']].drop_duplicates().sort_values(by='compositeLR_log10', ascending=False).reset_index(drop=True).reset_index().rename(columns=({'index': 'diseaseRank'})).drop(columns=['compositeLR_log10'])
    rank_df = rank_df.merge(gene_rank_df).merge(disease_rank_df)

    # Increment rank columns
    rank_df['geneRank'] += 1
    rank_df['diseaseRank'] += 1

    for gene in case_data['genes']:
         for d in case_data['genes'][gene]['diseases']:
            d.update({
                'geneRank': int(rank_df.loc[rank_df['omimId'] == d['omim'], 'geneRank'].reset_index(drop=True)[0]),
                'diseaseRank': int(rank_df.loc[rank_df['omimId'] == d['omim'], 'diseaseRank'].reset_index(drop=True)[0])
            })

    return case_data
