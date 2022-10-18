import pandas as pd
import os
import sys
import yaml
import json
import argparse

sys.path.append('..')
sys.path.append('.')

from config import *
from utils import *

# def run_lirical(case_id, clinphen_df, vcf, hpo_total, output_dir, tsv=True):
def run_lirical(case):
    
    # Filter ClinPhen tsv
    filtered_df = clinphen_df.head(hpo_total)

    # Isolate HPO IDs
    hpo_ids = filtered_df['HPO ID'].to_list()

    # Set temp and final filenames
    if tsv:
        temp_filename = os.path.join(config['project_root'], output_dir, case_id + '.tsv')
        output_filename = os.path.join(config['project_root'], output_dir, case_id + ".lirical.tsv")
    
    else:
        temp_filename = os.path.join(config['project_root'], output_dir, case_id + '.html')
        output_filename = os.path.join(config['project_root'], output_dir, case_id + ".lirical.html")
    
    try:
        
        # Write filtered df
        filtered_df.to_csv(os.path.join(config['project_root'], output_dir, case_id + "." + str(hpo_total)  + "TOTAL_filtered.clinphen.tsv"), sep='\t', index=False)

        # Build YAML
        yaml = build_yaml(template_path = os.path.join(config['project_root'], 'example', "template.yaml"), case_id = case_id, vcf = vcf, hpo_ids = hpo_ids, output_dir = output_dir, tsv=tsv)

        if not(os.path.exists(output_filename)):

            # Run LIRICAL
            os.system('cd ' + config['project_root'] + ' && java -jar LIRICAL.jar yaml -y ' + yaml)
            
            # Rename results files
            os.rename(temp_filename, output_filename)

        return output_filename
    
    except:
        print("Unable to run lirical, file reading error")
        sys.exit(1)




if __name__ == "__main__":

    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--hpo_list', '-hl', type=str, help='Path to the ordered list of HPO IDs.')
    parser.add_argument('--output', '-o', type=str, help='Path to a directory where the raw LIRICAL output data will be written.')
    parser.add_argument('--case_id', '-c', type=str, help='Case ID of the subject')
    parser.add_argument('--vcf', type=str, help='Path to case vcf')

    args = parser.parse_args()

    # Set variables
    case_id = args.case_id
    output_dir = args.output
    hpo_total = config['hpo_total_upper_bound']

    # Read in ClinPhen results
    clinphen_df = pd.read_csv(args.hpo_list, sep='\t', header=0)

    # Run html LIRICAL
    # run_lirical(case_id=case_id, clinphen_df=clinphen_df, hpo_total=hpo_total, output_dir=output_dir, tsv=False)
    run_lirical(case_id=case_id, clinphen_df=clinphen_df, vcf=args.vcf, hpo_total=len(clinphen_df.index), output_dir=output_dir, tsv=False)

    # Run tsv LIRICAL
    # run_lirical(case_id=case_id, clinphen_df=clinphen_df, hpo_total=hpo_total, output_dir=output_dir, tsv=True)
    run_lirical(case_id=case_id, clinphen_df=clinphen_df, vcf=args.vcf, hpo_total=len(clinphen_df.index), output_dir=output_dir, tsv=True)