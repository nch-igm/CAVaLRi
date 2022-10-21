import os
import pandas as pd

def get_diagnostic_df(case):

    config = case.cohort.config

    wes_diagnostic_df = pd.read_excel(os.path.join(case.cohort.root_path, config['diagnostic_source']))
    wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'] == 1]
    # wes_diagnostic_df = wes_diagnostic_df.loc[wes_diagnostic_df['Section'].notna()]
    wes_diagnostic_df = wes_diagnostic_df[['CoPath.Number', 'Disease.Gene.and.Transcript']]
    wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CoPath.Number':'case_id', 'Disease.Gene.and.Transcript':'gene'}))
    wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)
    return wes_diagnostic_df


def add_tp(case):

    # Get true positive data
    tp_df = get_diagnostic_df(case)
    
    # Intialize case hit list
    case.case_data['tpHits'] = []

    # Limit tp_df to only variants specific to this case
    tp_df = tp_df.loc[tp_df['case_id'] == case.case_id]
    
    # Check if the case has a diagnosis
    case.case_data['isTp'] = 1 if len(tp_df.index) != 0 else 0

    # Iterate through each disease
    for d in case.case_data['diseases']:
        
        # Check if the disease is related to a true positive gene
        if len(tp_df.loc[tp_df['gene'] == d['gene_data']['gene']].index) != 0:
            d['isHit'] = 1
            case.case_data['tpHits'].append(d['geneRank'])
        else:
            d['isHit'] = 0

    return case.case_data

