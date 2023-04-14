import sys
for path in ['.','../..']:
    sys.path.append(path)
import pickle
import argparse
import json
import os


def main(input, output):

    # Add source path
    root_path = os.path.abspath(os.path.dirname(__file__))
    root_path = root_path[:root_path.find('/src')]
    sys.path.append(root_path)

    # Read in cohort
    with open(input, 'rb') as f:
        case = pickle.load(f)

    # Process variants
    case.genotype.annotate_variants()
    case.genotype.filter_variants()
    case.genotype.read_variants()
    case.genotype.variants.to_csv('/Users/rsrxs003/projects/CAVaLRi_/variants.csv', index = False)

    # Build case data
    case.build_case_data()

    # Process phenotypes
    case.phenotype.read_phenotypes()

    case.phenotype.score_phenotypes()

    # Run LIRICAL
    # case.run_lirical()

    # Run case level mode of inheritance commands
    case.phenotype.calculate_phenotypeLR()
    case.genotype.get_scoring_annotations()
    case.genotype.get_pathogenic_variants()
    case.genotype.calculate_genoLR()
    case.calculate_moiLR()

    # # Run case level summary commands
    case.calculate_compositeLR()
    case.add_rankings()
    case.add_tp()
    case.build_case_summary()

    with open(output, 'wb') as f:
        pickle.dump(case, file = f)


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input', '-i', type=str, help='Case object to process')
    parser.add_argument('--output', '-o', type=str, help='Path to output case directory')
    args = parser.parse_args()

    main(args.input, args.output)
