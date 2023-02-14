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
    case.conda_bin = os.path.join(sys.exec_prefix, 'bin')
    case.genotype.annotate_variants()
    case.genotype.filter_variants()
    case.genotype.read_variants()

    # Process phenotypes
    case.phenotype.read_phenotypes()

    # Run LIRICAL
    case.run_lirical()

    # Run case level mode of inheritance commands
    case.build_case_data()
    case.phenotype.score_phenotypes()
    case.phenotype.calculate_phenotypeLR()
    case.calculate_genoLR()
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
