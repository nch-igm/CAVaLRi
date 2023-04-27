import math
import json
from scipy.stats import norm, skewnorm
import pickle


def get_diagnostic_probability(cLR, prior, pos_mu, pos_std, a, loc, scale):
    pos_prob = norm.pdf(cLR, pos_mu, pos_std)
    neg_prob = skewnorm.pdf(cLR, a, loc, scale)
    odds = (pos_prob/neg_prob) * prior
    probability = odds / (1 + odds)
    return probability

def calculate_compositeLR(case):

    config = case.cohort.config

    # For each OMIM ID, calculate the compositeLR
    for g in case.case_data['genes'].keys():

        for d, d_data in case.case_data['genes'][g].items():

            # Append genoLR
            d_data['genoLR'] = case.genotype.genotype_LR[g]

            # Append moiLR
            d_data['moiLR'] = case.moiLRs[d]
            
            # Append compositeLR
            d_data['compositeLR'] = d_data['phenoLR'] * config['phenoLR_scalar'] + d_data['genoLR'] * config['genoLR_scalar'] + d_data['moiLR'] * config['moiLR_scalar']

            # Calculate post-test probability
            posttest_probability = get_diagnostic_probability(
                    cLR = d_data['compositeLR'], 
                    prior = 0.006708, 
                    pos_mu = 10.47,
                    pos_std = 8.21,
                    a = 5.18,
                    loc = -17.47,
                    scale = 9.9
                )

            # Append post test probability
            d_data['postTestProbability'] =  posttest_probability

    # Return the result data frame
    return case.case_data
