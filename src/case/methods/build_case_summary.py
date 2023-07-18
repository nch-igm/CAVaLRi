import pandas as pd
import json

def build_case_summary(case):

    config = case.cohort.config

    # Initialize result data frame
    res_df = pd.DataFrame(columns=['geneRank', 'geneNcbiId', 'geneSymbol', 'postTestProbability', 'phenoLR', 'geneLR', 'moiLR', 'cavalriScore','hpoCount','variants'])

    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():
            if d != 'gene_data':
                print(g_data)
                res_df = pd.concat([
                    res_df, 
                    pd.DataFrame({
                        'geneRank': d_data['geneRank'],
                        'geneNcbiId': g,
                        'geneSymbol': g_data['gene_data']['symbol'],
                        'postTestProbability': d_data['postTestProbability'],
                        'phenoLR': d_data['phenoLR'] * config['phenoLR_scalar'],
                        'geneLR': g_data['gene_data']['genoLR'] * config['genoLR_scalar'],
                        'moiLR': d_data['moiLR'] * config['moiLR_scalar'],
                        'cavalriScore': d_data['score'],
                        'hpoCount': d_data['phenoCount'],
                        'variants': ';'.join([ v['hg38_position'] for v in g_data['gene_data']['variants'] ])
                    }, index = [0])
                ])

    # Get highest compositeLR for each disease
    compositLR_max_df = res_df[['geneNcbiId', 'postTestProbability']].groupby('geneNcbiId').max().reset_index().rename(columns=({'postTestProbability': 'max_postTestProbability'}))
    res_df = res_df.merge(compositLR_max_df)
    res_df = res_df.loc[res_df['postTestProbability'] == res_df['max_postTestProbability']].reset_index(drop=True).drop(columns=['max_postTestProbability']).drop_duplicates()

    return res_df.sort_values('geneRank')
