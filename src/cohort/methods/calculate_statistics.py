import os
import sys
import pandas as pd
import numpy as np
import scipy
from scipy import stats


def calculate_wilcoxon(summ_data, comp_data):

    # Limit to true positive cases
    summ_data = summ_data.loc[summ_data['isTp'] == 1]
    comp_data = comp_data.loc[comp_data['isTp'] == 1]
    
    # Intialize a new data frame to add ranks to
    rank_df = summ_data[['subjectId', 'hitRank']].rename(columns=({'hitRank': 'newRank'}))
    rank_df = rank_df.merge(comp_data[['subjectId', 'hitRank']].rename(columns=({'hitRank': 'compRank'})))
    
    # Remove all rows where there is a null ranking for either cohort
    rank_df = rank_df.loc[(rank_df['newRank'].notna()) & (rank_df['compRank'].notna())]
    
    # Calculate the one-sided Wilcoxon statistic
    data = scipy.stats.wilcoxon(rank_df['compRank'], rank_df['newRank'], alternative='greater')
    
    # Return the p-value
    return data.pvalue
    

def calculate_statistics(cohort):

    config = cohort.config

    try:
        print(f"Comparing results stored in {config['comparator']}")
    except:
        print(f'No comparator indicated')
        return 'No comparator'

    # Read in comparator data set
    summ_data = cohort.cohort_summary.copy()
    comp_data = pd.read_csv(os.path.join(config['comparator'], 'summary', 'cohort_summary.tsv'), sep='\t')
    assessed_cases = set(summ_data['subjectId']).intersection(set(comp_data['subjectId']))

    # Limit to shared cases
    summ_data = summ_data.loc[summ_data['subjectId'].isin(assessed_cases)]
    comp_data = comp_data.loc[comp_data['subjectId'].isin(assessed_cases)]

    # Calculate p-value of wilcoxon signed rank test 
    wilcoxon = calculate_wilcoxon(summ_data, comp_data)

    res_data = {
        'wilcoxon_p-value': wilcoxon#,
        # 'residuals': residuals
    }

    return res_data
