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
                            prior = 0.023376039559451566, 
                            pos_a = -2.5656189133285565,
                            pos_loc = 52.07742666463085,
                            pos_scale = 16.347160198789545,
                            neg_a = 30832529.791864797,
                            neg_loc = 6.863645946715076,
                            neg_scale = 11.328256872470757
                        )

                    if d_data['score'] > 61.735:
                        posttest_probability = 1

                if case.trio_status == 'DUO':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.014383989993746092, 
                            pos_a = 3.757650075261621,
                            pos_loc = 21.59155928397092,
                            pos_scale = 14.189580516852995,
                            neg_a = 1.9372801615618576,
                            neg_loc = 11.025804759643592,
                            neg_scale = 10.129620978366187
                        )

                    if d_data['score'] > 50.357:
                        posttest_probability = 0.5

                if case.trio_status == 'SINGLETON':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.017467248908296942, 
                            pos_a = -34917.468873275546,
                            pos_loc = 47.51852335253445,
                            pos_scale = 7.025608415541145,
                            neg_a = 38763923.22317009,
                            neg_loc = 13.38239685528061,
                            neg_scale = 7.58874516104736
                        )

                    if d_data['score'] > 47.517:
                        posttest_probability = 0.9

                # Append post test probability
                d_data['postTestProbability'] =  posttest_probability

    # Return the result data frame
    return case.case_data
