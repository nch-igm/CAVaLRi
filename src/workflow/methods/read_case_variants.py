import sys
sys.path.append('..')
sys.path.append('.')
import os
import pickle
from cavalri import *
import argparse

def main(input, case):

    # Read in cohort
    with open(input, 'rb') as f:
        x = pickle.load(f)

    # Find case
    for c in x.cases:
        if c.case_id == case:
            c.genotype.read_variants()
            return c


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--cohort', '-i', type=str, help='Cohort pickle')
    parser.add_argument('--case', '-c', type=str, help='Case being run')
    parser.add_argument('--output', '-o', type=str, help='Path to output case pickle')

    args = parser.parse_args()

    c = main(args.cohort, args.case)
    
    # Save the cohort
    with open(args.output, 'wb') as f:
        pickle.dump(c, file = f)
