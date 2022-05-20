import os
import sys
import pandas as pd
import json
import argparse

# Add package locations to sys.path
sys.path.append('..')
sys.path.append('.')

# import local config package
from config import *

def get_diagnostic_df():

    wes_diagnostic_df = pd.read_excel(os.path.join(config['project_root'], config['diagnostic_source']))
    wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'] == 1]
    # wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'].notna()]
    wes_diagnostic_df = wes_diagnostic_df[['CoPath.Number', 'Disease.Gene.and.Transcript']]

    for idx, row in wes_diagnostic_df.iterrows():
        gene = str(row['Disease.Gene.and.Transcript'])
        gene = gene[:gene.find('(')-1]
        wes_diagnostic_df.loc[idx, 'Disease.Gene.and.Transcript'] = gene

    wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CoPath.Number':'subjectId', 'Disease.Gene.and.Transcript':'gene'}))
    wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)
    return wes_diagnostic_df


def add_tp(case):

    # Get true positive data
    tp_df = get_diagnostic_df()
    print(tp_df)

    # Intialize case hit list
    case.case_data['tpHits'] = []

    # Limit tp_df to only variants specific to this case
    tp_df = tp_df.loc[tp_df['subjectId'] == case.case_id]
    
    # Check if the case has a diagnosis
    if len(tp_df.index) != 0:
        case.case_data['isTp'] = 1
    else:
        case.case_data['isTp'] = 0

    # Iterate through each gene
    for gene in case.case_data['genes']:

        # Check if the disease is related to a true positive gene
        if len(tp_df.loc[tp_df['gene'] == case.case_data['genes'][gene]['SYMBOL']].index) != 0:
            case.case_data['genes'][gene]['isHit'] = 1
            case.case_data['tpHits'].append(case.case_data['genes'][gene]['SYMBOL'])
        else:
            case.case_data['genes'][gene]['isHit'] = 0

    return case.case_data
