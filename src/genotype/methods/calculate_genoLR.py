import sys
import os
import pandas as pd
import json
import math
from scipy.stats import poisson
import time
fo = '/Users/rsrxs003/projects/CAVaLRi_/catch_some_output.txt'


def calculate_genoLR(genotype):
    """
    input:
        genotype -- CAVaLRi Genotype type object 
    output:
        dictionary keyed by gene symbol with a float value corresponding to the 
        log_10 scaled genotype likelihood ratio
    """

    # Read in the configuration settings of the Cohort object
    config = genotype.case.cohort.config

    # Intialize result dict
    res = {}

    for g in genotype.variants['GENE'].unique():
        
        g_df = genotype.variants[genotype.variants['GENE'] == g]

        # Calculate cumulative scores between all passing variants
        # score = g_df['conservation_score'].sum()
        score = len(g_df.index)
    
        # # Get maximum pathogenicity score
        # max_path_score = 0
        # for variant in d['gene_data']['variants']:
        #     if variant['pathScore'] > max_path_score:
        #         max_path_score = variant['pathScore']

        # # Calculate genotype likelihood ratio
        # if d['gene_data']['background_freq'] != 0.001 and d['gene_data']['disease_freq'] != 0.001:
        #     background_freq = 10**-5 + d['gene_data']['background_freq']
        #     disease_freq = d['gene_data']['disease_freq']
        #     # geneLR_log10 = math.log(poisson.pmf(disease_freq, max_path_score * disease_freq)/poisson.pmf(disease_freq, max_path_score * background_freq),10)
        #     geneLR_log10 = 0

        # else:
        #     geneLR_log10 = 0

        # Append to result
        res[g] = score

    return res