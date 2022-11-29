import sys
import os
import pandas as pd
import json
import math
from scipy.stats import poisson



def calculate_genoLR(case):

    config = case.cohort.config

    # Intialize result dict
    res = {'subjectId': case.case_id}

    # For each disease, access gene data and calculate genotype likelihood ratio

    for d in case.case_data['diseases']:
        
        # Get maximum pathogenicity score
        max_path_score = 0
        for variant in d['gene_data']['variants']:
            if variant['pathScore'] > max_path_score:
                max_path_score = variant['pathScore']

        # Calculate genotype likelihood ratio
        if d['gene_data']['background_freq'] != 0.001 and d['gene_data']['disease_freq'] != 0.001:
            background_freq = 10**-5 + d['gene_data']['background_freq']
            disease_freq = d['gene_data']['disease_freq']
            # geneLR_log10 = math.log(poisson.pmf(disease_freq, max_path_score * disease_freq)/poisson.pmf(disease_freq, max_path_score * background_freq),10)
            geneLR_log10 = 0

        else:
            geneLR_log10 = 0

        # Append to result
        res.update({
            # d['omimId']: geneLR_log10
            d['omimId']: d['gene_data']['scraped_geneLR']
        })

    return res