import os
import pandas as pd
import json
import time
fo = '/Users/rsrxs003/projects/CAVaLRi_/catch_some_output.txt'

def add_rankings(case):

    config = case.cohort.config

    # Intialize rank data frame
    rank_df = pd.DataFrame(columns=['omimId', 'geneSymbol', 'compositeLR'])
    
    # Populate the rank data frame
    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():

            # Create a dataframe to sort by compositeLR
            rank_df = pd.concat([
                rank_df,
                pd.DataFrame({
                'omimId': d,
                'geneSymbol': g,
                'compositeLR': d_data['compositeLR']
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
    
    # Assign values to the case.case_data attribute
    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():

            try:
                gene_rank = int(rank_df.loc[rank_df['omimId'] == d, 'geneRank'].reset_index(drop=True)[0])
            except:
                print(f'Checkpoint: 14 {g} {d} Calculating cLR {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
            
            disease_rank = int(rank_df.loc[rank_df['omimId'] == d, 'diseaseRank'].reset_index(drop=True)[0])

            d_data['geneRank'] = gene_rank # Step can be duplicated, but oh well
            d_data['diseaseRank'] = disease_rank

    # for d in case.case_data['diseases']:
    #     d.update({
    #         'geneRank': int(rank_df.loc[rank_df['omimId'] == d['omimId'], 'geneRank'].reset_index(drop=True)[0]),
    #         'diseaseRank': int(rank_df.loc[rank_df['omimId'] == d['omimId'], 'diseaseRank'].reset_index(drop=True)[0])
    #     })


    return case.case_data
