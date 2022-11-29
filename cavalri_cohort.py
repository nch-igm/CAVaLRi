from src.cohort import Cohort
from config import *
import argparse
import json
import shlex
import subprocess

def worker(cmd, err = False):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return err.decode() if err else out.decode()

def main(input_dir, output_dir):

    # Parse the input file
    output_dir = input_dir if not output_dir else output_dir
    cohort = Cohort(input_dir, output_dir, config)

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
    args = parser.parse_args()

    main(args.input_dir, args.output_dir)