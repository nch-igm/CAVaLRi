from src.cohort import Cohort
from src.cohort import Case
from config import *
import sys
import os
import pandas as pd
import json
import vcf
import subprocess
import shlex
import argparse
import pickle
import uuid
import yaml


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')


def update_config(temp_pickle_path, config_template_path, 
                    root_dir, config_output_path):
    """
    inputs:
        temp_pickle_path (str) - filepath to a cavalri.Case object.
        workflow_config_path(str): path to a local yaml file that informs a 
        cavalri.Workflow instance.
    outputs:
        Boolean, True if no errors were encountered, False otherwise.
    """

    # Read in the yaml
    with open(config_template_path,'r') as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
    
    # Add data
    cfg['root_dir'] = root_dir
    cfg['temp_pickle'] = temp_pickle_path
    cfg['temp_dir'] = os.path.dirname(temp_pickle_path)
    cfg['case'] = True

    # Overwrite yaml
    with open(config_output_path,'w') as f:
        cfg = yaml.dump(cfg, f)


def main(input, output_dir):

    # Parse the input file
    output_dir = os.path.dirname(input) if not output_dir else output_dir
    cohort = Cohort(os.path.dirname(input), output_dir, '', config)

    # Add case
    with open(input,'r') as d:
        inputs = json.load(d)
    case = os.path.basename(input)
    case = case[:case.find('.json')]

    # required_keys = ['phenotype','vcf','biological_sex','proband','mother_affected','father_affected']
    required_keys = ['phenotype','vcf','proband','pedigree']

    for rk in required_keys:
        if rk not in inputs.keys():
            print(f'Required key: {rk} was not provided in the CAVaLRi input file')
            sys.exit(1)

    # Read pedigree to get necessary inputs
    try:
        ped_df = pd.read_csv(inputs['pedigree'], header = None, sep = '\t')
        ped_df.columns = ['family','individual','father','mother','biological_sex','affected']
        proband_row = ped_df[ped_df['individual'] == inputs['proband']].reset_index(drop=True)
        inputs['biological_sex'] = 'M' if proband_row.loc[0,'biological_sex'] == 1 else 'F'
        inputs['mother'] = proband_row.loc[0,'mother']
        inputs['father'] = proband_row.loc[0,'father']
        for parent in ['mother','father']:
            if inputs[parent] == 0:
                inputs[parent] = 'Unavailable'
                inputs[f'{parent}_affected'] = 0
            else:
                parent_row = ped_df[ped_df['individual'] == inputs[parent]].reset_index(drop=True)
                inputs[f'{parent}_affected'] = 1 if parent_row.loc[0,'affected'] == 2 else 0

    except Exception as e:
        print(f"{type(e)}: {e}")
        print("Ped file not formatted correctly, see https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format")
        sys.exit(1)


    # # Check for parents
    # for parent in ['mother','father']:
    #     if parent in inputs.keys():
    #         if inputs[parent] == '':
    #             inputs[parent] = 'Unavailable'
    #     else:
    #         inputs[parent] = 'Unavailable'
    

    # Intialize case object
    cs = Case(
        cohort = cohort, 
        case_id = case,
        phenotype_path = inputs['phenotype'],
        genotype_path = inputs['vcf'],
        biological_sex = inputs['biological_sex'],
        proband = inputs['proband'],
        mother = inputs['mother'],
        father = inputs['father'],
        mother_affected = inputs['mother_affected'],
        father_affected = inputs['father_affected']
    )

    # Set up a temporary directory
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

    # Write temp file path to workflow config.yaml
    workflow_path = os.path.join(cohort.root_path, 'src/workflow')
    workflow_config_path = os.path.join(workflow_path, 'config.yaml')
    config_output_path = os.path.join(cs.temp_dir, 'config.yaml')
    update_config(case_pickle_path, workflow_config_path, cohort.root_path, config_output_path)
    
    # Run snakemake pipeline
    cmd = f"cd {workflow_path} && {os.path.join(conda_bin, 'snakemake')} --cores {config['cores']} --configfile {config_output_path}"
    p = worker(cmd)

    # Check to see if the pipeline ran successfully
    if not os.path.exists(full_pickle_path):
        print(p)
        sys.exit(1)

    # Load result
    with open(full_pickle_path, 'rb') as f:
        cs = pickle.load(f)

    # Add scored phenotypes to case data
    cs.case_data['phenotypes'] = cs.phenotype.phenotypes
    cs.case_data['genes'] = {int(k):v for k,v in cs.case_data['genes'].items()}

    # Write output files
    with open(os.path.join(output_dir, f'{cs.case_id}.cavalri.json'), 'w') as f:
        json.dump(cs.case_data, f, indent = 4)
    cs.case_summary.to_csv(os.path.join(output_dir, f'{cs.case_id}.cavalri.summary.csv'), index = False)

    # Remove temporary directory
    worker(f'rm -Rf {temp_folder}')

    # Get versions of dependencies
    # hpo_version = worker(f"grep #date: {os.path.join(cs.cohort.root_path,  config['lirical_data_path'], 'phenotype.hpoa')}")
    # config['HPOA version'] = hpo_version.split(' ')[-1].strip()

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