import sys
sys.path.append('../..')
import pickle
import argparse

from config import *

def main(input, output):

    # Read in cohort
    with open(input, 'rb') as f:
        case = pickle.load(f)

    # Process variants
    case.genotype.read_variants()
    case.genotype.annotate_variants()
    # case.genotype.filter_variants()

    # Process phenotypes
    case.phenotype.read_phenotypes()

    # Run LIRICAL
    case.run_lirical()

    # Run case level mode of inheritance commands
    # case.build_case_data()
    # case.calculate_moiLR()

    # # Run case level phenotype commands
    # case.calculate_genotypeLR()
    # case.phenotype.read_phenotypes()
    # case.score_phenotypes()
    # case.calculate_phenotypeLR()

    # # Run case level summary commands
    # case.calculate_compositeLR()
    # case.add_rankings()
    # case.add_tp()
    # case.build_case_summary()

    with open(output, 'wb') as f:
        pickle.dump(case, file = f)


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input', '-i', type=str, help='Case object to process')
    parser.add_argument('--output', '-o', type=str, help='Path to output case directory')
    args = parser.parse_args()

    main(args.input, args.output)
