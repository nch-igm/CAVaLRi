import sys
import math
import pandas as pd



def calculate_phenotypeLR(case):

    config = case.cohort.config

    for gene, gene_diseases in case.case_data['genes'].items():
        for d, d_data in gene_diseases.items():

            try:

                # Initialize a data frame to capture composite likelihood ratios for each iteration
                # run_df = pd.DataFrame(columns = ['hpo_total', 'phenoLR'])
                hpo_totals, composite_phenoLRs = [], []

                # Get the necessary phenotype rank and LR data
                pheno_df = pd.DataFrame(d_data['phenotype_scores']).T
                pheno_rank_df = pd.DataFrame(case.phenotype.phenotypes).T
                pheno_rank_df = pheno_rank_df.merge(pheno_df, left_index = True, right_index = True)[['rank', 'LR']]
                
                for hpo_total in range(config['hpo_lower_bound'], config['hpo_upper_bound'] + 1):
                    
                    composite_phenoLR = pheno_rank_df.loc[pheno_rank_df['rank'] <= hpo_total]['LR'].product()

                    # Add HPO term count, compositeLR and post-test probability to the run data frame
                    # run_df = run_df.append({
                    #     'hpo_total': hpo_total,
                    #     'phenoLR': composite_phenoLR
                    # }, ignore_index = True)

                    hpo_totals.append(hpo_total)
                    composite_phenoLRs.append(composite_phenoLR)
                
                run_df = pd.DataFrame({'hpo_total':hpo_totals, 'phenoLR':composite_phenoLRs})

                # Get maximum values of compositeLR and total phenotypes considered from the run data frame
                max_phenoLR = run_df['phenoLR'].max()
                                
                # Get positions where the max phenoLR is occuring
                phenoCount = [str(i+config['hpo_lower_bound']) for i, j in enumerate(run_df['phenoLR']) if j == max_phenoLR]
                phenoCount = ",".join(phenoCount)

                # Append OMIM ID, PhenoLR and associated HPO set length to result
                gene_diseases[d]['phenoLR_log10'] = round(math.log(max_phenoLR, 10), 3)
                gene_diseases[d]['phenoCount'] = phenoCount
                            
            except Exception as e:
                # print(f"An exception: {e} occured at {time.strftime('%H:%M:%S', time.localtime())}: {d}", file = open(fo, 'a'))
                gene_diseases[d]['phenoLR_log10'] = 0
                gene_diseases[d]['phenoCount'] = None

    # Return the result data frame
    return case.case_data
