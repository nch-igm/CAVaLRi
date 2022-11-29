import pandas as pd

def is_tp(row, hits):
    return 1 if row['geneRank'] in hits else 0

def compile_cohort_variants(cohort):
    dfs = []
    for case in cohort.cases:
        df = case.case_summary.copy()
        df['case'] = case.case_id
        df['result'] = df.apply(is_tp, hits = case.case_data['tpHits'], axis = 1)
        dfs.append(df)
    return pd.concat(dfs).sort_values('compositeLR', ascending = False)
    