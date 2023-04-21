import sys
import pickle
import argparse
import json
import os


def main(input, output):

    # Add source path
    root_path = os.path.abspath(os.path.dirname(__file__))
    root_path = root_path[:root_path.find('/src')]
    sys.path.append(root_path)

    # Read in pheno case
    with open(input, 'rb') as f:
        pheno_case = pickle.load(f)

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

    # Run aggregate methods
    pheno_case.calculate_compositeLR()
    pheno_case.add_tp()
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
