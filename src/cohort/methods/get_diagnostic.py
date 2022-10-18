import pandas as pd
import os

def get_diagnostic_df(cohort):

    wes_diagnostic_df = pd.read_excel(os.path.join(cohort.root_path, cohort.config['diagnostic_source']))
    wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'] == 1]
    wes_diagnostic_df = wes_diagnostic_df[['CoPath.Number', 'Disease.Gene.and.Transcript']]

    for idx, row in wes_diagnostic_df.iterrows():
        gene = str(row['Disease.Gene.and.Transcript'])
        if gene.find('(') != -1:
            gene = gene[:gene.find('(')-1]
        wes_diagnostic_df.loc[idx, 'Disease.Gene.and.Transcript'] = gene

    wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CoPath.Number':'subjectId', 'Disease.Gene.and.Transcript':'gene'}))
    wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)
    return wes_diagnostic_df


def get_extended_diagnostic_df(cohort):

    wes_diagnostic_df = pd.read_excel(os.path.join(cohort.root_path, cohort.config['diagnostic_source']))
    wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'] == 1]
    # wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'].notna()]
    # wes_diagnostic_df = wes_diagnostic_df[['CoPath.Number', 'Disease.Gene.and.Transcript']]

    for idx, row in wes_diagnostic_df.iterrows():
        gene = str(row['Disease.Gene.and.Transcript'])
        gene = gene[:gene.find('(')-1]
        wes_diagnostic_df.loc[idx, 'Disease.Gene.and.Transcript'] = gene

    wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CoPath.Number':'subjectId', 'Disease.Gene.and.Transcript':'gene'}))
    wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)
    return wes_diagnostic_df