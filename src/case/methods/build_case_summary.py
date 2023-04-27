import pandas as pd

def build_case_summary(case):

    config = case.cohort.config

    # Initialize result data frame
    # res_df = pd.DataFrame(columns=['geneRank', 'geneNcbiId', 'geneSymbol', 'postTestProbability', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR','hpoCount'])
    res_df = pd.DataFrame(columns=['geneRank', 'geneSymbol', 'postTestProbability', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR','hpoCount'])

    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():
            res_df = pd.concat([
                res_df, 
                pd.DataFrame({
                    'geneRank': d_data['geneRank'],
                    # 'geneNcbiId': d['gene_data']['geneNcbiId'],
                    'geneSymbol': g,
                    'postTestProbability': d_data['postTestProbability'],
                    'phenoLR': d_data['phenoLR'] * config['phenoLR_scalar'],
                    'geneLR': d_data['genoLR'] * config['genoLR_scalar'],
                    'moiLR': d_data['moiLR'] * config['moiLR_scalar'],
                    'compositeLR': d_data['compositeLR'],
                    'hpoCount': d_data['phenoCount']
                }, index = [0])
            ])

    # Get highest compositeLR for each disease
    compositLR_max_df = res_df[['geneSymbol', 'compositeLR']].groupby('geneSymbol').max().reset_index().rename(columns=({'compositeLR': 'max_compositeLR'}))
    res_df = res_df.merge(compositLR_max_df)
    res_df = res_df.loc[res_df['compositeLR'] == res_df['max_compositeLR']].reset_index(drop=True).drop(columns=['max_compositeLR']).drop_duplicates()

    return res_df.sort_values('geneRank')
