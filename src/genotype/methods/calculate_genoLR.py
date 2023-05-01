import sys
import os
import pandas as pd
import json
import math
from scipy.stats import poisson


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
    root_path = genotype.case.cohort.root_path

    # Intialize result dict
    res = {}

    # Read in mode of inheritance
    moi_df = pd.read_csv(os.path.join(root_path, config['moi_db']))
    disease_count = len(moi_df.index)
    dominant_codes = set(['AD','XLD'])
    recessive_codes = set(['AR','XLR'])

    pathogenic_df = genotype.pathogenic_variants.copy()

    for g in list(genotype.case.case_data['genes'].keys()):
        
        # Get alt allele count
        g_df = pathogenic_df[pathogenic_df['GENE'] == g]
        def get_gt(row):
            gt = json.loads(row['proband'].replace("'",'\"'))['GT']
            if gt in ['0|1','1/1','1|1'] or (gt in ['0/1','1/0'] and genotype.case.biological_sex == 'M' and row['CHROM'] == 'X'):
                return 2
            elif gt in ['0/1','1/0']:
                return 1
            else:
                return 0
        
        g_df['var_count'] = g_df.apply(get_gt, axis = 1)
        var_count = g_df['var_count'].sum()

        # Determine mode of inheritance
        diseases = [f'OMIM:{d}' for d in genotype.case.case_data['genes'][g].keys()]
        mois = list(moi_df[moi_df['omimId'].isin(diseases)]['moi'].unique())
        mois = ';'.join(mois)
        mois = set(mois.split(';'))

        # Compare variant count to mode of inheritance
        if len(dominant_codes & mois) >= 1 and var_count >= 1:
            score = g_df['score'].max()
        elif len(recessive_codes & mois) >= 1 and var_count >= 2:
            score = g_df['score'].mean()
        else:
            score = 0

        # Multiply score by the length of possible diseases
        score = score * disease_count

        # Square the value if a ClinVar pathogenic is present
        clinvar_count = len(g_df[g_df['clinvar_path_sig']].index)
        if clinvar_count > 0:
            score = score ** 2

        res[g] = math.log(score, 10)

    return res