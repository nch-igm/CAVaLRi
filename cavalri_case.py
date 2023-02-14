from src.cohort import Cohort
from src.cohort import Case
from config import *
import sys
import os
import pandas as pd
import json
import subprocess
import shlex
import argparse
import pickle
import uuid

def worker(cmd, err = False):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return err.decode() if err else out.decode()


def validate_case(inputs):

    path_variables = ['phenotype','vcf']
    for pv in path_variables:
        if not os.path.exists(inputs[pv]):
            return pv
    return ''

def main(input, output_dir):

    # Parse the input file
    output_dir = os.path.dirname(input) if not output_dir else output_dir
    cohort = Cohort(os.path.dirname(input), output_dir, '', config)

    # Add case
    with open(input,'r') as d:
        inputs = json.load(d)
    case = os.path.basename(input)
    case = case[:case.find('.json')]

    # Check for parents
    for parent in ['mother','father']:
        if parent in inputs.keys():
            if inputs[parent] == '':
                inputs[parent] = 'Unavailable'
        else:
            inputs[parent] = 'Unavailable'

    validated_input = validate_case(inputs)
    if validated_input != '':
        print(f"Validation failed: {case}, {inputs[validated_input]} not a valid path")

    else:
        cs = Case(
            cohort = cohort, 
            case_id = case,
            phenotype_path = inputs['phenotype'],
            genotype_path = inputs['vcf'],
            biological_sex = inputs['biological_sex'],
            proband = inputs['proband'],
            mother = inputs['mother'],
            father = inputs['father']
        )

        # Set up the temporary directory
        temp_folder = os.path.abspath(os.path.join(output_dir, str(uuid.uuid4())))
        cs.temp_dir = temp_folder
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)

        # cohort.add_case(cs)

        # Pickle the case
        case_pickle_path = os.path.join(temp_folder, f'{cs.case_id}.pickle')
        with open(case_pickle_path, 'wb') as f:
            pickle.dump(cs, file = f)

        # Run case
        full_pickle_path = os.path.join(temp_folder, f'{cs.case_id}.full.pickle')
        script_path = os.path.join(cohort.root_path, 'src/workflow/scripts/run_case.py')
        conda_bin = os.path.join(sys.exec_prefix, 'bin')
        p = worker(f"{os.path.join(conda_bin, 'python')} {script_path} --input {case_pickle_path} -o {full_pickle_path}", err=True)

        if not os.path.exists(full_pickle_path):
            print(p)
            sys.exit(1)

        # Load result
        with open(full_pickle_path, 'rb') as f:
            cs = pickle.load(f)
        with open(os.path.join(output_dir, f'{cs.case_id}.cavalri.json'), 'w') as f:
            json.dump(cs.case_data, f, indent = 4)
        cs.case_summary.to_csv(os.path.join(output_dir, f'{cs.case_id}.cavalri.summary.csv'), index = False)

        # Remove temporary directory
        worker(f'rm -Rf {temp_folder}')

        # Get LIRICAL and HPOA version
        lirical_version = worker(f"java -jar {os.path.join(cs.cohort.root_path,  config['lirical_executable'])} --version")
        hpo_version = worker(f"grep #date: {os.path.join(cs.cohort.root_path,  config['lirical_data_path'], 'phenotype.hpoa')}")
        config['LIRICAL version'] = lirical_version.strip()
        config['HPOA version'] = hpo_version.split(' ')[-1].strip()

        # Create run log
        log_path = os.path.join(output_dir,'config.json')
        with open(log_path, 'w') as f:
            json.dump(config, f, indent=4)

        

if __name__ == '__main__':

    # Intialize parser
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--input', '-i', type=str, 
        help='Path to CAVaLRi subject input file')
    parser.add_argument('--output_dir', '-o', type=str, 
        help='Directory where CAVaLRi output files are written')
    
    args = parser.parse_args()
    main(args.input, args.output_dir)