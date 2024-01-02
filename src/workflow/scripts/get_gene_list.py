import sys
import pickle
import subprocess
import argparse
import json
import os
import time


def worker(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    return out.decode() if out else err.decode()

def vcf_is_empty(vcf_path):
    cmd = f"gunzip -c {vcf_path} | grep -vc '^#'"
    count = worker(cmd)
    return int(count) == 0

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
    case.genotype.split_variants()
    is_empty = vcf_is_empty(case.genotype.short_variant_path)
    if not is_empty:
        case.genotype.normalize_variants()
        case.genotype.annotate_variants()
        case.genotype.filter_variants()
        case.genotype.read_filtered_variants()
    case.genotype.score_cnvs()
    case.genotype.score_pathogenicity()
    case.genotype.filter_pathogenic_variants(is_empty = is_empty)

    # Build case data
    case.build_case_data()

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
