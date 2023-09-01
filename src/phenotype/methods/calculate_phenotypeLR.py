import sys
import math
import pandas as pd
import pickle
import json
import re


def calculate_phenotypeLR(case):

    config = case.cohort.config

    for gene, gene_diseases in case.case_data['genes'].items():
        for d, d_data in gene_diseases.items():

            if not re.search('gene_data', d):

                try:

                    hpo_totals, composite_phenoLRs = [], []

                    # Get the necessary phenotype rank and LR data
                    pheno_df = pd.DataFrame(d_data['phenotype_scores']).T
                    pheno_df = pheno_df[pheno_df['common_ancestor'].notna()]
                    pheno_rank_df = pd.DataFrame(case.phenotype.phenotypes).T
                    pheno_rank_df = pheno_rank_df.merge(pheno_df, left_index = True, right_index = True)[['rank', 'LR']]

                    # Transform to log space
                    pheno_rank_df['log(LR)'] = pheno_rank_df['LR'].apply(math.log, args = (10,))
                    
                    # for hpo_total in range(config['hpo_lower_bound'], config['hpo_upper_bound'] + 1):
                    for hpo_total in range(config['hpo_lower_bound'], len(pheno_rank_df.index) + 1):
                        
                        composite_phenoLR = pheno_rank_df.loc[pheno_rank_df['rank'] <= hpo_total]['log(LR)'].sum()
                        hpo_totals.append(hpo_total)
                        composite_phenoLRs.append(composite_phenoLR)
                    
                    run_df = pd.DataFrame({'hpo_total':hpo_totals, 'phenoLR':composite_phenoLRs})

                    # Get maximum values of compositeLR and total phenotypes considered from the run data frame
                    max_phenoLR = run_df['phenoLR'].max()
                                    
                    # Get positions where the max phenoLR is occuring
                    phenoCount = [str(i+config['hpo_lower_bound']) for i, j in enumerate(run_df['phenoLR']) if j == max_phenoLR]
                    phenoCount = ",".join(phenoCount)

                    # Append OMIM ID, PhenoLR and associated HPO set length to result
                    gene_diseases[d]['phenoLR'] = round(max_phenoLR, 3)
                    gene_diseases[d]['phenoCount'] = phenoCount
                                
                except Exception as e:
                    # print(f"An exception: {e} occured at {time.strftime('%H:%M:%S', time.localtime())}: {d}", file = open(fo, 'a'))
                    gene_diseases[d]['phenoLR'] = 0
                    gene_diseases[d]['phenoCount'] = 0

    # Return the result data frame
    return case.case_data
