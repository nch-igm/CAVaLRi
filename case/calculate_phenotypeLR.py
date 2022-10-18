import sys
import math
import pandas as pd

sys.path.append('..')
# from hpo_walk import ontology
# from hpo_walk.annotations import get_phenotype_disease_gene_df
from config import *

def calculate_phenotypeLR(case):

    for gene in case.case_data['genes']:
        for d in case.case_data['genes'][gene]['diseases']:

            # Initialize a data frame to capture composite likelihood ratios for each iteration
            run_df = pd.DataFrame(columns = ['hpo_total', 'phenoLR'])

            # Get the necessary phenotype rank and LR data
            pheno_df = pd.DataFrame(index = d['phenotype_scores'].keys(), data = d['phenotype_scores'].values())
            pheno_rank_df = pd.DataFrame(index = case.phenotype.phenotypes.keys(), data = {'rank': case.phenotype.phenotypes.values()})
            pheno_rank_df = pheno_rank_df.merge(pheno_df, left_index = True, right_index = True)[['rank', 'LR']]


            for hpo_total in range(config['hpo_total_lower_bound'], config['hpo_total_upper_bound'] + 1):
                
                composite_phenoLR = pheno_rank_df.loc[pheno_rank_df['rank'] <= hpo_total]['LR'].product()

                # Add HPO term count, compositeLR and post-test probability to the run data frame
                run_df = run_df.append({
                    'hpo_total': hpo_total,
                    'phenoLR': composite_phenoLR
                }, ignore_index = True)

            # Get maximum values of compositeLR and total phenotypes considered from the run data frame
            max_phenoLR = run_df['phenoLR'].max()
            
            # Get positions where the max phenoLR is occuring
            phenoCount = [str(i+config['hpo_total_lower_bound']) for i, j in enumerate(run_df['phenoLR']) if j == max_phenoLR]
            phenoCount = ",".join(phenoCount)

            # Append OMIM ID, PhenoLR and associated HPO set length to result
            d['phenoLR_log10'] = round(math.log(max_phenoLR, 10), 3)
            d['phenoCount'] = phenoCount

    # Return the result data frame
    return case.case_data
