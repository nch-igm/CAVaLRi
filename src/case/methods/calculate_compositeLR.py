import math
import json
from scipy.stats import norm, skewnorm
import pickle
import re


def get_diagnostic_probability(score, prior, pos_a, pos_loc, pos_scale, neg_a, neg_loc, neg_scale):
    try:
        pos_prob = skewnorm.pdf(score, pos_a, pos_loc, pos_scale)
        neg_prob = skewnorm.pdf(score, neg_a, neg_loc, neg_scale)
        odds = (pos_prob/neg_prob) * prior
        probability = odds / (1 + odds)
        return probability
    except:
        return 0.000001

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
                            prior = 0.028203723986856513,
                            pos_a = -3.196236533168645,
                            pos_loc = 52.90631610920187,
                            pos_scale = 16.42447592324006,
                            neg_a = 88522229.93833573,
                            neg_loc = 6.86373631034402,
                            neg_scale = 11.79271054053001
                        )

                    if d_data['score'] > 61.197:
                        posttest_probability = 1

                if case.trio_status == 'DUO':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.014383989993746092,
                            pos_a = 3.3995746521502754,
                            pos_loc = 21.47500644802009,
                            pos_scale = 13.637011889082167,
                            neg_a = 1.5048478711331637,
                            neg_loc = 10.776897326128825,
                            neg_scale = 9.940154478558602
                        )

                    if d_data['score'] > 50.5:
                        posttest_probability = 0.7

                if case.trio_status == 'SINGLETON':
                    posttest_probability = get_diagnostic_probability(
                            score = d_data['score'], 
                            prior = 0.017467248908296942,
                            pos_a = 0.2745585486276849,
                            pos_loc = 39.87043047425281,
                            pos_scale = 3.869089531065029,
                            neg_a = 28459576.362230994,
                            neg_loc = 13.43428827583272,
                            neg_scale = 6.836695812458304
                        )

                    if d_data['score'] > 46.296:
                        posttest_probability = 0.9

                # Append post test probability
                d_data['postTestProbability'] =  posttest_probability

    # Return the result data frame
    return case.case_data
