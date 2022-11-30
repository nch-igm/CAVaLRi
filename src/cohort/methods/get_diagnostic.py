import pandas as pd
import os

def get_diagnostic_df(cohort):

    wes_diagnostic_df = pd.read_csv(case.cohort.diagnostic_data)
    wes_diagnostic_df = wes_diagnostic_df[['CASE', 'DIAGNOSTIC_GENE']]
    wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CASE':'case_id', 'DIAGNOSTIC_GENE':'gene'}))

    # for idx, row in wes_diagnostic_df.iterrows():
    #     gene = str(row['Disease.Gene.and.Transcript'])
    #     if gene.find('(') != -1:
    #         gene = gene[:gene.find('(')-1]
    #     wes_diagnostic_df.loc[idx, 'Disease.Gene.and.Transcript'] = gene

    wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)
    return wes_diagnostic_df


def get_extended_diagnostic_df(cohort):

    wes_diagnostic_df = pd.read_csv(case.cohort.diagnostic_data)
    wes_diagnostic_df = wes_diagnostic_df[['CASE', 'DIAGNOSTIC_GENE']]
    wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CASE':'case_id', 'DIAGNOSTIC_GENE':'gene'}))

    wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)
    return wes_diagnostic_df