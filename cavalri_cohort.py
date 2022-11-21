from src.cohort import Cohort
from config import *
import argparse

def main(input_dir, output_dir):

    # Parse the input file
    output_dir = input_dir if not output_dir else output_dir
    cohort = Cohort(input_dir, output_dir, config)

    # Run CAVaLRi workflow for each case in the cohort
    cohort.run()

if __name__ == '__main__':

    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_dir', '-i', type=str, help='Directory where CAVaLRi subject input files are stored')
    parser.add_argument('--output_dir', '-o', type=str, help='Directory where CAVaLRi output files are written')
    args = parser.parse_args()

    main(args.input_dir, args.output_dir)