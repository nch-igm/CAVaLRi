import math
import json
from scipy.stats import norm, skewnorm
import pickle
import re


def get_diagnostic_probability(score, prior, pos_a, pos_loc, pos_scale, neg_a, neg_loc, neg_scale):
    pos_prob = skewnorm.pdf(score, pos_a, pos_loc, pos_scale)
    neg_prob = skewnorm.pdf(score, neg_a, neg_loc, neg_scale)
    odds = (pos_prob/neg_prob) * prior
    probability = odds / (1 + odds)
    return probability

def calculate_compositeLR(case):

    config = case.cohort.config

    # For each OMIM ID, calculate the compositeLR
    for g, g_data in case.case_data['genes'].items():

        for d, d_data in g_data.items():

            if not re.search('gene_data', d):
                
                # Append compositeLR
                d_data['score'] = d_data['phenoLR'] * config['phenoLR_scalar'] + g_data['gene_data']['genoLR'] * config['genoLR_scalar'] + d_data['moiLR'] * config['moiLR_scalar']
                
                # Calculate post-test probability based on trio status
                if case.trio_status == 'TRIO':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.02375, 
                            pos_a = -4.3068,
                            pos_loc = 58.9998,
                            pos_scale = 20.9783,
                            neg_a = 56669107.1977,
                            neg_loc = 0.7469,
                            neg_scale = 13.7409
                        )

                if case.trio_status == 'DUO':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.01436, 
                            pos_a = .0000093,
                            pos_loc = 34.0425,
                            pos_scale = 10.0607,
                            neg_a = 0.6955,
                            neg_loc = 11.6084,
                            neg_scale = 11.5504
                        )

                if case.trio_status == 'SINGLETON':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.0165, 
                            pos_a = -265905.7711,
                            pos_loc = 46.52007,
                            pos_scale = 7.1597,
                            neg_a = 13618624.1039,
                            neg_loc = 12.8368,
                            neg_scale = 9.3558
                        )

                # Append post test probability
                d_data['postTestProbability'] =  posttest_probability

    # Return the result data frame
    return case.case_data
