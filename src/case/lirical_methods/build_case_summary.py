import pandas as pd

def build_case_summary(case):

    config = case.cohort.config

    # Initialize result data frame
    # res_df = pd.DataFrame(columns=['geneRank', 'geneNcbiId', 'geneSymbol', 'postTestProbability', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR','hpoCount'])
    res_df = pd.DataFrame(columns=['geneRank', 'geneSymbol', 'postTestProbability', 'phenoLR', 'geneLR', 'moiLR', 'compositeLR','hpoCount'])

    for d in case.case_data['diseases']:
        res_df = pd.concat([
            res_df, 
            pd.DataFrame({
                'geneRank': d['geneRank'],
                # 'geneNcbiId': d['gene_data']['geneNcbiId'],
                'geneSymbol': d['gene_data']['gene'],
                'postTestProbability': d['postTestProbability'],
                'phenoLR': d['phenoLR'] * config['phenoLR_scalar'],
                'geneLR': d['genoLR'] * config['genoLR_scalar'],
                'moiLR': d['moiLR'] * config['moiLR_scalar'],
                'compositeLR': d['compositeLR'],
                'hpoCount': d['hpoCount']
            }, index = [0])
        ])

    # Get highest compositeLR for each disease
    compositLR_max_df = res_df[['geneSymbol', 'compositeLR']].groupby('geneSymbol').max().reset_index().rename(columns=({'compositeLR': 'max_compositeLR'}))
    res_df = res_df.merge(compositLR_max_df)
    res_df = res_df.loc[res_df['compositeLR'] == res_df['max_compositeLR']].reset_index(drop=True).drop(columns=['max_compositeLR']).drop_duplicates()

    return res_df.sort_values('geneRank')
