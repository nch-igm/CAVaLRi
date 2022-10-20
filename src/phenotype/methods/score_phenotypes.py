import os
import pandas as pd


def score_phenotypes(phenotype):

    config = phenotype.case.cohort.config

    # Read in clinphen tsv to get HPO ID rank
    all_pheno_scores = pd.read_csv(os.path.join(phenotype.case.cohort.root_path, config['pheno_score_source']))
    all_pheno_scores = all_pheno_scores[['HPO ID', config['pheno_score']]]
    pheno_case_df = phenotype.pheno_df.merge(all_pheno_scores, on = 'HPO ID')
    pheno_case_df = pheno_case_df.sort_values(config['pheno_score'], ascending = False).reset_index(drop=True)
    pheno_case_df = pheno_case_df.reset_index().rename(columns={'index': 'phenotypeRank'})
    pheno_case_df['phenotypeRank'] +=1
    pheno_case_df = pheno_case_df[['HPO ID', 'phenotypeRank']]

    # Create a disease level attribute to store HPO Rankings
    phenotype.case.case_data.update({
        'hpoRankings': []
    })

    # Get the rank for each HPO term
    for hpo in phenotype.case.case_data['diseases'][0]['pheno_data']:
        # print(hpo)
        pheno_rank = pheno_case_df.loc[pheno_case_df['HPO ID'] == hpo['hpoId']].reset_index(drop=True).loc[0, 'phenotypeRank']
        phenotype.case.case_data['hpoRankings'].append({
            'hpoId': hpo['hpoId'],
            'hpoRank': int(pheno_rank)
        })
