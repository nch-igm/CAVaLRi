import os
import pandas as pd
import json


def add_rankings(case):

    config = case.cohort.config

    # Intialize rank data frame
    rank_df = pd.DataFrame(columns=['omimId', 'geneSymbol', 'compositeLR'])

    for d in case.case_data['diseases']:
        rank_df = pd.concat([
            rank_df,
            pd.DataFrame({
            'omimId': d['omimId'],
            'geneSymbol': d['gene_data']['gene'],
            'compositeLR': d['compositeLR']
        }, index = [0])
    ])

    # Build disease and gene ranks
    gene_rank_df = rank_df[['geneSymbol', 'compositeLR']].groupby('geneSymbol').max().reset_index().sort_values(by='compositeLR', ascending=False)
    gene_rank_df = gene_rank_df.reset_index(drop=True).reset_index().rename(columns=({'index': 'geneRank'})).drop(columns=['compositeLR'])
    disease_rank_df = rank_df[['omimId', 'compositeLR']].drop_duplicates().sort_values(by='compositeLR', ascending=False).reset_index(drop=True).reset_index().rename(columns=({'index': 'diseaseRank'})).drop(columns=['compositeLR'])
    rank_df = rank_df.merge(gene_rank_df).merge(disease_rank_df)

    # Increment rank columns
    rank_df['geneRank'] += 1
    rank_df['diseaseRank'] += 1

    for d in case.case_data['diseases']:
        d.update({
            'geneRank': int(rank_df.loc[rank_df['omimId'] == d['omimId'], 'geneRank'].reset_index(drop=True)[0]),
            'diseaseRank': int(rank_df.loc[rank_df['omimId'] == d['omimId'], 'diseaseRank'].reset_index(drop=True)[0])
        })


    return case.case_data
