import sys
sys.path.append('..')
sys.path.append('.')
import os
import pickle
from cavalri import *
import argparse

def main(case_pickle_dir, cohort):

    # Read in cohort
    with open(cohort, 'rb') as f:
        x = pickle.load(f)

    # Iterate through case variant pickles
    for entry in os.scandir(case_pickle_dir):
        if entry.name.endswith('.pickle'):
            
            # Read case pickle and append to cohort
            with open(os.path.join(case_pickle_dir, entry.name), 'rb') as f:
                c = pickle.load(f)

            # Add to cohort
            case = entry.name[:entry.name.find('.variants.pickle')]
            for xc in x.cases:
                if xc.case_id == case:
                    xc.genotype.variants = c.genotype.variants

    # Run cohort annotation commands
    x.aggregate_variants()
    print('Variants aggregated')
    x.send_unique_variants()
    print('Variants sent to S3')
    x.run_snpeff()
    print('SNPeff complete')    
    x.run_annotation_notebook()
    print('EMR notebook executed')
    x.retrieve_annotated_variants()
    print('Annotated variants retrieved')

    # Save the cohort
    with open(cohort, 'wb') as f:
        pickle.dump(x, file = f)

    x.score_variants()
    print('Variants scored')

    return x


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--case_pickle_dir', '-i', type=str, help='Case directory containing variant pickle')
    parser.add_argument('--cohort', '-c', type=str, help='Cohort pickle')
    parser.add_argument('--output', '-o', type=str, help='Path to output case directory')

    args = parser.parse_args()

    x = main(args.case_pickle_dir, args.cohort)

    # Create output file
    os.system('touch {}'.format(args.output))

    # Save the cohort
    with open(args.cohort, 'wb') as f:
        pickle.dump(x, file = f)

