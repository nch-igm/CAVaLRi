import os
import pandas as pd
import pickle

def get_diagnostic_df(case):

    config = case.cohort.config
    
    if case.cohort.diagnostic_data != '':
        wes_diagnostic_df = pd.read_csv(case.cohort.diagnostic_data)
        wes_diagnostic_df = wes_diagnostic_df[['CASE', 'DIAGNOSTIC_GENE']]
        wes_diagnostic_df = wes_diagnostic_df.rename(columns=({'CASE':'case_id', 'DIAGNOSTIC_GENE':'gene'}))
        wes_diagnostic_df = wes_diagnostic_df.drop_duplicates().reset_index(drop=True)

    else:
        wes_diagnostic_df = pd.DataFrame(columns = ['case_id', 'gene'])

    return wes_diagnostic_df


def add_tp(case):

    with open('/Users/rsrxs003/projects/CAVaLRi_/case.pickle','wb') as f:
        pickle.dump(case, f)

    # Get true positive data
    tp_df = get_diagnostic_df(case)
    
    # Intialize case hit list
    case.case_data['tpHits'] = []

    # Limit tp_df to only variants specific to this case
    tp_df = tp_df.loc[tp_df['case_id'] == case.case_id]
    
    # Check if the case has a diagnosis
    case.case_data['isTp'] = 1 if len(tp_df.index) != 0 else 0

    # Iterate through each disease
    for g, g_data in case.case_data['genes'].items():
        for d, d_data in g_data.items():
        
            # Check if the disease is related to a true positive gene
            if len(tp_df.loc[tp_df['gene'] == g].index) != 0:
                d_data['isHit'] = 1
                case.case_data['tpHits'].append(d_data['geneRank'])
            else:
                d_data['isHit'] = 0

    return case.case_data

