import sys
for path in ['.','../..']:
    sys.path.append(path)
import pickle
import argparse
import json
import os
import time
fo = '/Users/rsrxs003/projects/CAVaLRi_/catch_some_output.txt'


def main(input, output):

    # Add source path
    root_path = os.path.abspath(os.path.dirname(__file__))
    root_path = root_path[:root_path.find('/src')]
    sys.path.append(root_path)
    print(f'Checkpoint: 4 {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))


    # Read in cohort
    with open(input, 'rb') as f:
        case = pickle.load(f)

    # Process variants
    case.genotype.annotate_variants()
    case.genotype.filter_variants()
    case.genotype.read_variants()
    print(f'Checkpoint: 5 {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
    case.genotype.variants.to_csv('/Users/rsrxs003/projects/CAVaLRi_/variants.csv', index = False)

    # Build case data
    case.build_case_data()
    print(f'Checkpoint: 6 {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))


    # Process phenotypes
    case.phenotype.read_phenotypes()
    print(f'Checkpoint: 7 {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))

    case.phenotype.score_phenotypes()
    print(f'Checkpoint: 10 {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
    with open('/Users/rsrxs003/projects/CAVaLRi_/catch_some_phenos.json', 'w') as f:
        json.dump(case.phenotype.scored_phenotypes, f, indent=4)

    # Run LIRICAL
    # case.run_lirical()

    # Run case level mode of inheritance commands
    case.phenotype.calculate_phenotypeLR()
    print(f'Checkpoint: 11 {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
    case.genotype.get_scoring_annotations()
    case.genotype.get_pathogenic_variants()
    case.genotype.calculate_genoLR()
    print(f'Checkpoint: 12 Calculating MOI {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
    case.calculate_moiLR()

    # # Run case level summary commands
    print(f'Checkpoint: 13 Calculating cLR {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
    case.calculate_compositeLR()
    case.add_rankings()
    case.add_tp()

    pickle_path = '/Users/rsrxs003/projects/CAVaLRi_/case.pickle'
    with open(pickle_path, 'wb') as f:
        pickle.dump(case, f)
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
