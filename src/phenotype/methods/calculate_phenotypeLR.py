import json
import pandas as pd


def calculate_phenotypeLR(case):

    config = case.cohort.config

    # Get HPO Ranks
    hpo_rank_df = pd.DataFrame(data = case.case_data['hpoRankings']).sort_values(by='hpoRank').reset_index(drop=True)

    # Intialize result dict
    res = {
        'subjectId': case.case_id
    }

    # For each OMIM ID, find the maximium phenotype likelihood ratio

    for d in case.case_data['diseases']:

        # Initialize OMIM specific data frame to capture each HPO term list phenoLR
        run_df = pd.DataFrame(columns = ['hpo_total', 'phenoLR'])

        # Iterate through each HPO term lists length and calculate compositeLR and post-test probability,
        # then take the max of each.  We also need to grab the length of the list when the max was achieved.

        for hpo_total in range(config['hpo_lower_bound'], config['hpo_upper_bound'] + 1):
        # for hpo_total in range(len(hpo_rank_df.index), len(hpo_rank_df.index) + 1):
            
            # Sum together hpoLRs, grouped by omimId
            omim_pheno_df = pd.DataFrame(data = d['pheno_data'])
            omim_pheno_df = hpo_rank_df.merge(omim_pheno_df)
            # print(omim_pheno_df)
            phenoLR = omim_pheno_df.loc[0:hpo_total-1]['hpoLR'].sum()

            # Add HPO term count, compositeLR and post-test probability to the run data frame
            run_df = pd.concat([run_df, pd.DataFrame({
                'hpo_total': hpo_total,
                'phenoLR': phenoLR
            }, index = [0])])

        # Get maximum values of compositeLR and postTestProbability from the run data frame
        # print(d['gene_data']['gene'])
        # print(omim_pheno_df)
        # print(run_df)
        max_phenoLR = run_df['phenoLR'].max()
        # print(max_phenoLR)
        # Get positions where the max phenoLR is occuring
        hpoCount = [str(i+config['hpo_lower_bound']) for i, j in enumerate(run_df['phenoLR']) if j == max_phenoLR]
        hpoCount = ",".join(hpoCount)

        # Append OMIM ID, PhenoLR and associated HPO set length toresult
        res.update({
            d['omimId']: {
                'phenoLR': round(max_phenoLR, 3),
                'hpoCount': hpoCount
            }
        })

    # Return the result data frame
    return res