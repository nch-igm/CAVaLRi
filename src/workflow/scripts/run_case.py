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

    # Read in cohort
    with open(input, 'rb') as f:
        case = pickle.load(f)

    # Process variants
    case.validate_inputs()
    case.genotype.annotate_variants()
    case.genotype.filter_variants()
    case.genotype.read_filtered_variants()
    case.genotype.score_pathogenicity()
    case.genotype.filter_pathogenic_variants()
    case.genotype.get_pathogenic_variants()

    # Build case data
    case.build_case_data()

    # Process phenotypes
    case.phenotype.read_phenotypes()
    case.phenotype.score_phenotypes()
    case.phenotype.calculate_phenotypeLR()

    # Process genotype likelihoods
    case.genotype.calculate_genoLR()

    # Process mode of inheritence likelihoods
    case.calculate_moiLR()

    # Run aggregate methods
    case.calculate_compositeLR()
    case.add_tp()
    case.add_rankings()
    case.build_case_summary()

    # Save the processed case object
    with open(output, 'wb') as f:
        pickle.dump(case, file = f)


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input', '-i', type=str, help='Case object to process')
    parser.add_argument('--output', '-o', type=str, help='Path to output case directory')
    args = parser.parse_args()

    main(args.input, args.output)
