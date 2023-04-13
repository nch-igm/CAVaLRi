import math
import json
from scipy.stats import norm, skewnorm
import pickle
import time
fo = '/Users/rsrxs003/projects/CAVaLRi_/catch_some_output.txt'


def get_diagnostic_probability(cLR, prior, pos_mu, pos_std, a, loc, scale):
    pos_prob = norm.pdf(cLR, pos_mu, pos_std)
    neg_prob = skewnorm.pdf(cLR, a, loc, scale)
    odds = (pos_prob/neg_prob) * prior
    probability = odds / (1 + odds)
    return probability

def calculate_compositeLR(case):

    config = case.cohort.config

    # For each OMIM ID, calculate the compositeLR
    # for g in case.case_data['genes'].keys():
    # case.genotype.pathogenic_variants.to_csv('/Users/rsrxs003/projects/CAVaLRi_/pathogenic_variants.csv', index = False)

    for g in case.case_data['genes'].keys():

        # with open('/Users/rsrxs003/projects/CAVaLRi_/case_data_.json','w') as f:
        #     json.dump(case.case_data, fp = f, indent = 4)

        for d, d_data in case.case_data['genes'][g].items():

            # Append phenoLR
            # case.case_data['genes'][g][d].update({
            #     'phenoLR': case.phenotype.phenoLRs[d['omimId']]['phenoLR'],
            #     'hpoCount': case.phenotype.phenoLRs[d['omimId']]['hpoCount']
            # })

            # Append genoLR
            d_data['genoLR'] = case.genotype.genotype_LR[g]

            # Append moiLR
            d_data['moiLR'] = case.moiLRs[d]
            
            # Append compositeLR
            try:
                d_data['compositeLR'] = d_data['phenoLR_log10'] * config['phenoLR_scalar'] + d_data['genoLR'] * config['genoLR_scalar'] + d_data['moiLR'] * config['moiLR_scalar']
            except:
                print(f'ERROR {g} {d} {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))

            

            # Calculate post-test probability
            # numerator = d['pretest_probability'] * 10**d['compositeLR']
            # denominator = (1 - d['pretest_probability']) + d['pretest_probability'] * 10**d['compositeLR']
            # posttest_probability = numerator / denominator
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
