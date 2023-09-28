from config import *
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
import CAVaLRi as cv
import pandas as pd
import json
import subprocess
import argparse
import yaml


def worker(cmd):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')



def update_config(temp_pickle_path: str, config_template_path: str, 
                    root_dir: str, config_output_path: str):
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


def main(input_dir, output_dir, diagnostic_data):

    # Validate input and output directories
    output_dir = input_dir if not output_dir else output_dir
    for dir in [input_dir,output_dir]:
        try:
            if not os.path.exists(dir):
                os.mkdir(dir)
        except:
            print(f'{dir} is not a valid directory')
            sys.exit(1)
    
    cohort = cv.Cohort(os.path.dirname(input), output_dir, config)

    cohort.make_temp_dir()

    required_keys = ['phenotype','vcf','proband','pedigree']

    # Add cases
    for entry in os.scandir(input_dir):
         
        # Intialize remove case flag
        rm_case = 0

        # Read in input file
        with open(os.path.join(input_dir, entry.name),'r') as d:
            inputs = json.load(d)
        case = os.path.basename(input)
        case = case[:case.find('.json')]

        # Check required keys
        for rk in required_keys:
            if rk not in inputs.keys():
                print(f'Required key: {rk} was not provided in the CAVaLRi input file')
                print(f'Removing {case}')
                rm_case = 1
                continue
        if rm_case == 1:
            continue


        # Read pedigree to get necessary inputs
        try:
            ped_df = pd.read_csv(inputs['pedigree'], header = None, sep = '\t')
            ped_df.columns = ['family','individual','father','mother','biological_sex','affected']
            ped_df = ped_df.astype({x:str for x in ['individual','father','mother']})
            proband_row = ped_df[ped_df['individual'] == inputs['proband']].reset_index(drop=True)
            inputs['biological_sex'] = 'M' if proband_row.loc[0,'biological_sex'] == 1 else 'F'
            inputs['mother'] = proband_row.loc[0,'mother']
            inputs['father'] = proband_row.loc[0,'father']
            for parent in ['mother','father']:
                if inputs[parent] == '0':
                    inputs[parent] = 'Unavailable'
                    inputs[f'{parent}_affected'] = 0
                else:
                    parent_row = ped_df[ped_df['individual'] == inputs[parent]].reset_index(drop=True)
                    inputs[f'{parent}_affected'] = 1 if parent_row.loc[0,'affected'] == 2 else 0

        except Exception as e:
            print(f"{type(e)}: {e}")
            print("Ped file not formatted correctly, see https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format")
            print(f'Removing {case}')
            continue
        

        # Intialize case object
        cs = cv.Case(
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

        cohort.add_case(cs)
    
    # Validate diagnostic data
    # if diagnostic_data:
    #     diagnostic_data = os.path.abspath(diagnostic_data)
    #     diagnostic_df = pd.read_csv(diagnostic_data)
    #     if len(set(diagnostic_df.columns).intersection(set(['CASE','DIAGNOSTIC_GENE']))) != 2:
    #         print(f"'CASE' and 'DIAGNOSTIC_GENE' columns were not found in {diagnostic_data}")


    # Run CAVaLRi workflow for each case in the cohort
    cohort.run()

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