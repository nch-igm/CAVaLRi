from config import *
import sys
import os
import pickle
import pandas as pd
import json
import subprocess
import argparse
import yaml
sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
import CAVaLRi as cv

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


def main(input_dir, output_dir):

    # Validate input and output directories
    output_dir = input_dir if not output_dir else output_dir
    case_output_dir = os.path.join(output_dir, 'cases')

    for dir in [input_dir,output_dir,case_output_dir]:
        try:
            if not os.path.exists(dir):
                os.mkdir(dir)
        except:
            print(f'{dir} is not a valid directory')
            sys.exit(1)
    
    cohort = cv.Cohort(input_dir, output_dir, config)

    cohort.make_temp_dir()

    required_keys = ['phenotype','vcf','proband','pedigree']

    # Add cases
    for entry in [ e for e in os.scandir(input_dir) if e.name.endswith('.json')]:
                 
        # Intialize remove case flag
        rm_case = 0

        # Read in input file
        with open(os.path.join(input_dir, entry.name),'r') as d:
            inputs = json.load(d)
        case = os.path.basename(entry.name)
        case = case[:case.find('.json')]

        # Check required keys
        for rk in required_keys:
            if rk not in inputs.keys():
                print(f'Required key: {rk} was not provided in the CAVaLRi input file')
                print(f'Removing {case}')
                rm_case = 1
                break
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

        # Pickle case
        case_pickle_path = os.path.join(cs.cohort.temp_dir, f'{cs.case_id}.pickle')
        with open(case_pickle_path, 'wb') as f:
            pickle.dump(cs, file = f)

        # Set finished pickle path
        cs.full_pickle_path = os.path.join(cs.cohort.temp_dir, f'{cs.case_id}.full.pickle')

        cohort.add_case(cs)

    
    # Validate diagnostic data
    # if diagnostic_data:
    #     diagnostic_data = os.path.abspath(diagnostic_data)
    #     diagnostic_df = pd.read_csv(diagnostic_data)
    #     if len(set(diagnostic_df.columns).intersection(set(['CASE','DIAGNOSTIC_GENE']))) != 2:
    #         print(f"'CASE' and 'DIAGNOSTIC_GENE' columns were not found in {diagnostic_data}")


    # Run cohort
    conda_bin = os.path.join(sys.exec_prefix, 'bin')

    # Write temp file path to workflow config.yaml
    workflow_path = os.path.join(cohort.root_path, 'src/workflow')
    workflow_config_path = os.path.join(workflow_path, 'config.yaml')
    config_output_path = os.path.join(cohort.temp_dir, 'config.yaml')
    update_config(case_pickle_path, workflow_config_path, cohort.root_path, config_output_path)

    print(f'Running CAVaLRi for {len(cohort.cases.keys())} cases: ({", ".join(cohort.cases.keys())})')
 
    # Run snakemake pipeline
    cmd = f"cd {workflow_path} && {os.path.join(conda_bin, 'snakemake')} --cores {config['cores']} --configfile {config_output_path}"
    p = worker(cmd)

    # Check to see if the pipeline ran successfully
    failed_cases = []
    successful_cases = []
    for k,v in cohort.cases.items():
        if not os.path.exists(v.full_pickle_path):
            failed_cases.append(k)
        else:
            successful_cases.append(k)
    
    if len(failed_cases) != 0:
        print(f'The following cases failed to execute: {";".join(failed_cases)}')
        for fc in failed_cases:
            cohort.remove_case(fc)

    # Create output
    cohort_summary = []

    for k,v in cohort.cases.items():

        # Load result
        with open(v.full_pickle_path, 'rb') as f:
            c = pickle.load(f)

        # Add scored phenotypes to case data
        c.case_data['phenotypes'] = c.phenotype.phenotypes
        c.case_data['genes'] = {int(k):v_ for k,v_ in c.case_data['genes'].items()}

        # Write output files
        with open(os.path.join(case_output_dir, f'{k}.cavalri.json'), 'w') as f:
            json.dump(c.case_data, f, indent = 4)
        c.case_summary.to_csv(os.path.join(case_output_dir, f'{k}.cavalri.summary.csv'), index = False)
        c.case_summary['case'] = k
        cohort_summary.append(c.case_summary)

    print(f'Writing output to {output_dir}')

    summary_path = os.path.join(output_dir, 'gene_summary.csv')
    pd.concat(cohort_summary).sort_values('postTestProbability', ascending=False).to_csv(summary_path, index = False)

    with open(os.path.join(output_dir, 'config.json'), 'w') as f:
        json.dump(config, f, indent=4)

    # Remove temporary directory
    cohort.remove_temp_dir()


if __name__ == '__main__':

    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_dir', '-i', type=str, help='Directory where CAVaLRi subject input files are stored')
    parser.add_argument('--output_dir', '-o', type=str, help='Directory where CAVaLRi output files are written')
    args = parser.parse_args()

    main(args.input_dir, args.output_dir)