import sys
import pickle
import argparse
import pandas as pd
import json
import os
import re


def main(input, output):

    # Add source path
    root_path = os.path.abspath(os.path.dirname(__file__))
    root_path = root_path[:root_path.find('/src')]
    sys.path.append(root_path)

    # Read in pheno case
    with open(input, 'rb') as f:
        pheno_case = pickle.load(f)

    config = pheno_case.cohort.config

    # Read in geno case
    geno_case_path = os.path.join(os.path.dirname(input), f'{pheno_case.case_id}.geno.pickle')
    with open(geno_case_path, 'rb') as f:
        geno_case = pickle.load(f)

    # Read in moi case
    moi_case_path = os.path.join(os.path.dirname(input), f'{pheno_case.case_id}.moi.pickle')
    with open(moi_case_path, 'rb') as f:
        moi_case = pickle.load(f)
    
    pheno_case.genotype = geno_case.genotype
    pheno_case.moiLRs = moi_case.moiLRs

    # Get diseases that lack mendelien inheritance annotations
    moi_df = pheno_case.cohort.moi.copy()
    no_moi_diseases = list(moi_df[moi_df['moi'].isna()]['omimId'].unique())

    # Limit cases to only those that have a phenotype likelihood ratio
    new_case_data = {'genes':{g:{} for g in pheno_case.case_data['genes'].keys()}}
    for g, g_data in new_case_data['genes'].items():

        for d, d_data in pheno_case.case_data['genes'][g].items():
            if 'phenoCount' in d_data.keys():
                if d_data['phenoCount'] != 0 and f'OMIM:{d}' not in no_moi_diseases:
                    g_data[d] = d_data
                    g_data[d]['moiLR'] = pheno_case.moiLRs[d]
            if re.search('gene_data', d):
                g_data[d] = d_data
                g_data[d]['genoLR'] = pheno_case.genotype.genotype_LR[g]
    
    pheno_case.case_data = {'genes':{k:v for k,v in new_case_data['genes'].items() if len(v.keys()) > 1 and pheno_case.genotype.genotype_LR[k] != 0}}

    # Run aggregate methods
    pheno_case.calculate_compositeLR()
    # pheno_case.add_tp()
    pheno_case.add_rankings()
    pheno_case.build_case_summary()

    # Save the processed case object
    with open(output, 'wb') as f:
        pickle.dump(pheno_case, file = f)


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input', '-i', type=str, help='Case object to process')
    parser.add_argument('--output', '-o', type=str, help='Path to output case directory')
    args = parser.parse_args()

    main(args.input, args.output)
