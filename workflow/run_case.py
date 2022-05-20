import sys
sys.path.append('..')
sys.path.append('.')
import pickle
from cavalri import *
import argparse

def main(input, case):

    # Read in cohort
    with open(input, 'rb') as f:
        x = pickle.load(f)

    # Get case of interest
    for c in x.cases:
        if c.case_id == case:
            case = c

    # Run case level genotype commands
    c.genotype.read_variants()
    c.genotype.annotate_variants()
    c.genotype.filter_variants()

    # Run case level mode of inheritance commands
    c.build_case_data()
    c.calculate_moiLR()

    # Run case level phenotype commands
    c.calculate_genotypeLR()
    c.phenotype.read_phenotypes()
    c.score_phenotypes()
    c.calculate_phenotypeLR()

    # Run case level summary commands
    c.calculate_compositeLR()
    c.add_rankings()
    c.add_tp()
    c.build_case_summary()

    return x


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--cohort', '-i', type=str, help='Cohort pickle')
    parser.add_argument('--case', '-c', type=str, help='Case being run')
    parser.add_argument('--output', '-o', type=str, help='Path to output case directory')

    args = parser.parse_args()

    x = main(args.cohort, args.case)

    with open(args.output, 'wb') as f:
        pickle.dump(x, file = f)
