import math
import json

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
        def get_lirical_phenoLR(data):
            sum = 0
            for pheno in data['pheno_data']:
                sum += pheno['hpoLR']
            return sum

        d.update({
            'compositeLR': get_lirical_phenoLR(d) * config['phenoLR_scalar'] + d['gene_data']['scraped_geneLR'] * config['genoLR_scalar']
        })

        # Calculate post-test probability
        numerator = d['pretest_probability'] * 10**d['compositeLR']
        denominator = (1 - d['pretest_probability']) + d['pretest_probability'] * 10**d['compositeLR']
        posttest_probability = numerator / denominator

        # Append post test probability
        d.update({'postTestProbability': posttest_probability})

    # Return the result data frame
    return case.case_data
