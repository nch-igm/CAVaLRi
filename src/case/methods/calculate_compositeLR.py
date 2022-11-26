import math
import json
from scipy.stats import norm, skewnorm

def get_diagnostic_probability(cLR, prior, pos_mu, pos_std, a, loc, scale):
    pos_prob = norm.pdf(cLR, pos_mu, pos_std)
    neg_prob = skewnorm.pdf(cLR, a, loc, scale)
    odds = (pos_prob/neg_prob) * prior
    probability = odds / (1 + odds)
    return probability

def calculate_compositeLR(case):

    config = case.cohort.config

    # For each OMIM ID, calculate the compositeLR
    for d in case.case_data['diseases']:

        # Append phenoLR
        d.update({
            'phenoLR': case.phenotype.phenoLRs[d['omimId']]['phenoLR'],
            'hpoCount': case.phenotype.phenoLRs[d['omimId']]['hpoCount']
        })

        # Append genoLR
        d.update({
            'genoLR': case.genoLRs[d['omimId']]
        })

        # Append moiLR
        d.update({
            'moiLR': case.moiLRs[d['omimId']]
        })

        # Append compositeLR
        d.update({
            'compositeLR': d['phenoLR'] * config['phenoLR_scalar'] + d['gene_data']['scraped_geneLR'] * config['genoLR_scalar'] + d['moiLR'] * config['moiLR_scalar']
        })

        # Calculate post-test probability
        # numerator = d['pretest_probability'] * 10**d['compositeLR']
        # denominator = (1 - d['pretest_probability']) + d['pretest_probability'] * 10**d['compositeLR']
        # posttest_probability = numerator / denominator
        posttest_probability = get_diagnostic_probability(
                cLR = d['compositeLR'], 
                prior = 0.006708, 
                pos_mu = 10.47,
                pos_std = 8.21,
                a = 5.18,
                loc = -17.47,
                scale = 9.9
            )

        # Append post test probability
        d.update({'postTestProbability': posttest_probability})

    # Return the result data frame
    return case.case_data
