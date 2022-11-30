from src.cohort import Cohort
from config import *
import sys
import argparse
import pandas as pd
import json
import shlex
import subprocess

def worker(cmd, err = False):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return err.decode() if err else out.decode()

def main(input_dir, output_dir, diagnostic_data):

    # Validate input and output directories
    output_dir = input_dir if not output_dir else output_dir
    for dir in [input_dir,output_dir]:
        if not os.path.isdir(dir):
            print(f'{dir} is not a valid directory')
            sys.exit(1)
    
    # Validate diagnostic data
    if diagnostic_data:
        diagnostic_data = os.path.abspath(diagnostic_data)
        diagnostic_df = pd.read_csv(diagnostic_data)
        if len(set(diagnostic_df.columns).intersection(set(['CASE','DIAGNOSTIC_GENE']))) != 2:
            print(f"'CASE' and 'DIAGNOSTIC_GENE' columns were not found in {diagnostic_data}")

    # Build the cohort
    cohort = Cohort(input_dir, output_dir, diagnostic_data, config)

    # Run CAVaLRi workflow for each case in the cohort
    cohort.run()

    # Get LIRICAL and HPOA version
    lirical_version = worker(f"java -jar {os.path.join(cohort.root_path,  config['lirical_executable'])} --version")
    hpo_version = worker(f"grep #date: {os.path.join(cohort.root_path,  config['lirical_data_path'], 'phenotype.hpoa')}")
    config['LIRICAL version'] = lirical_version.strip()
    config['HPOA version'] = hpo_version.split(' ')[-1].strip()

    # Create run log
    log_path = os.path.join(output_dir, 'summary', 'config.json')
    with open(log_path, 'w') as f:
        json.dump(config, f, indent=4)

if __name__ == '__main__':

    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_dir', '-i', type=str, help='Directory where CAVaLRi subject input files are stored')
    parser.add_argument('--output_dir', '-o', type=str, help='Directory where CAVaLRi output files are written')
    parser.add_argument('--diagnostic_data', '-d', type=str, help='File containing two columns, CASE and DIAGNOSTIC_GENE indicating diagnostic genes for cases provided in the cohort')
    args = parser.parse_args()

    main(args.input_dir, args.output_dir, args.diagnostic_data)