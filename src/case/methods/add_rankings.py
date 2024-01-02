import os
import re
import pandas as pd
import json

def add_rankings(case):

    config = case.cohort.config

    # Intialize rank data frame
    rank_df = pd.DataFrame(columns=['omimId', 'geneId', 'postTestProbability'])
    
    # Populate the rank data frame
    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():
            if not re.search('gene_data', d):

                # Create a dataframe to sort by postTestProbability
                rank_df = pd.concat([
                    rank_df,
                    pd.DataFrame({
                    'omimId': d,
                    'geneId': g,
                    'postTestProbability': d_data['postTestProbability']
                    }, index = [0])
                ])

    # Build disease and gene ranks
    gene_rank_df = rank_df[['geneId', 'postTestProbability']].groupby('geneId').max().reset_index().sort_values(by='postTestProbability', ascending=False)
    gene_rank_df = gene_rank_df.reset_index(drop=True).reset_index().rename(columns=({'index': 'geneRank'})).drop(columns=['postTestProbability'])
    disease_rank_df = rank_df[['omimId', 'postTestProbability']].drop_duplicates().sort_values(by='postTestProbability', ascending=False).reset_index(drop=True).reset_index().rename(columns=({'index': 'diseaseRank'})).drop(columns=['postTestProbability'])
    rank_df = rank_df.sort_values(by='postTestProbability', ascending=False).reset_index(drop=True).merge(gene_rank_df).merge(disease_rank_df)

    # Increment rank columns
    rank_df['geneRank'] += 1
    rank_df['diseaseRank'] += 1
    
    # Assign values to the case.case_data attribute
    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():
            if not re.search('gene_data', d):

                gene_rank = int(rank_df.loc[rank_df['geneId'] == g, 'geneRank'].reset_index(drop=True)[0])            
                disease_rank = int(rank_df.loc[rank_df['omimId'] == d, 'diseaseRank'].reset_index(drop=True)[0])

                d_data['geneRank'] = gene_rank # Step can be duplicated, but that's okay
                d_data['diseaseRank'] = disease_rank

    # return case.case_data
