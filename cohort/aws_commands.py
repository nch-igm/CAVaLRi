import os
import sys
from tkinter.tix import Tree
import pandas as pd
import shlex
import subprocess
import json
import time
from datetime import date
import re
sys.path.append('../..')
sys.path.append('..')
from config import *


def send_unique_variants(cohort):

    # Set unique variant data path
    unfiltered_variants_path = os.path.join(config['project_root'], config['temp_data'], config['s3_csv_input'])
    
    # Write the data frame
    cohort.unique_variants.to_csv(unfiltered_variants_path, index = False)

    # Send to s3
    os.system('aws s3 cp {} {}'.format(unfiltered_variants_path, os.path.join(config['s3_bucket'], config['s3_csv_input'])))


def retrieve_annotated_variants(cohort):

    # Make local directory where csvs will be saved
    annotated_path = os.path.join(config['project_root'], config['temp_data'], 'annotated_variants')
    if not os.path.exists(annotated_path):
        os.mkdir(annotated_path)
    else:
        os.system('rm -r {}'.format(annotated_path))

    # Read results from S3
    os.system('aws s3 cp --quiet --recursive {} {}'.format(
        os.path.join(config['s3_bucket'], config['s3_output']),
        annotated_path
    ))

    # Organize output, join together all of the csv partitions
    # for entry in os.scandir(annotated_path):
    #     if entry.name.endswith('csv'):
    #         df = pd.read_csv(os.path.join(annotated_path, entry.name), dtype=str)
    #         res_df = df if 'res_df' not in locals() else res_df.append(df)
    res_df = pd.read_parquet(annotated_path)
    
    # Change data type of numeric columns
    data_types = {
        'total_samples': float,
        'alt_sample_fraction': float,
        'gnomad_ex_faf95_popmax': float,
        'yprob': float,
        'dbnsfp_phylop17way_primate': float,
        'dbnsfp_phylop17way_primate_rankscore': float,
        'dbnsfp_phastcons17way_primate_rankscore': float,
        'dbnsfp_primateai_rankscore': float,
        'dbnsfp_primateai_score': float,
        'dbnsfp_phastcons17way_primate': float,
        'SURF': float,
        'SURF_Stability': float
    }

    for k, v in data_types.items():
        try:
            res_df[k] = res_df[k].astype(v)
        except:
            pass

    return res_df.reset_index(drop=True)


def run_annotation_notebook(cohort):

    # Empty out the output directory
    os.system("aws s3 rm --quiet --recursive {}".format(os.path.join(config['s3_bucket'], config['s3_output'])))

    # Run EMR notebook
    command = """
    aws emr --region %s \
        start-notebook-execution \
        --editor-id %s \
        --notebook-params '{"input":"%s", "output":"%s", "varhouse_anno_path": "%s"}' \
        --relative-path %s \
        --notebook-execution-name my-execution \
        --execution-engine '{"Id" : "%s"}' \
        --service-role EMR_Notebooks_DefaultRole \
    """ % (config['aws_region'],
           config['emr_editor-id'],
           os.path.join(config['s3_bucket'], config['s3_csv_input']),
           os.path.join(config['s3_bucket'], config['s3_output']),
        #    'CAVaLRi_SNPeff_{}'.format(date.today().isoformat()),
           config['snpeff_parquet'],
           config['emr_notebook_name'],
           config['emr_execution-engine'])
    cmd = shlex.split(command)
    p = subprocess.Popen(cmd,  stdout=subprocess.PIPE)
    out, err = p.communicate()

    # Parse out execution ID
    response = json.loads(out)
    execution_id = response['NotebookExecutionId']

    # Check on the EMR notebook execution
    command = """
    aws emr --region %s \
    describe-notebook-execution --notebook-execution-id %s
    """ % (config['aws_region'], execution_id)
    cmd = shlex.split(command)
    while True:
        p = subprocess.Popen(cmd,  stdout=subprocess.PIPE)
        out, err = p.communicate()
        response = json.loads(out)['NotebookExecution']
        if response['Status'] == 'FINISHED':
            break
        if response['Status'] == 'FAILED':
            print(response)
            break
        time.sleep(10)


def run_snpeff(cohort):

    # Indicate template vcf
    template_path = os.path.join(config['project_root'], config['vcf_template'])

    # Indicate output file name
    new_vcf = os.path.join(config['project_root'], config['temp_cohort_vcf'])
    with open(new_vcf, 'w') as fo:

        # Copy headers
        with open(template_path, 'r') as f:
            for line in f:
                if re.search('^#', line):
                    print(line.strip(), file = fo)

        # Write each diagnostic variant
        for idx, row in cohort.unique_variants.iterrows():
            print('chr' + row['CHROM'], row['POS'], '.', row['REF'], row['ALT'], '.', '.', '.', '.', '.', '.', '.', file = fo, sep='\t')

    # Send the vcf to s3
    os.system('aws s3 cp {} {}'.format(new_vcf, os.path.join(config['s3_bucket'], config['s3_vcf_input'])))

    # Set up source directory
    os.system('ssh -i {} hadoop@{} "mkdir {}"'.format(config['igm_pem'], config['emr_cluster_ip'], config['emr_root']))

    # Create bash script to run on cluster
    bash_template = os.path.join(config['project_root'], config['varhouse_bash'])
    new_bash = os.path.join(config['project_root'], config['temp_data'], os.path.basename(config['varhouse_bash']))
    with open(new_bash, 'w') as fo:
        with open(bash_template, 'r') as f:

            # Define parameter update values
            map = {
                    '{output_definition}': os.path.join(config['emr_root'], os.path.basename(config['varhouse_output_definition'])),
                    '{analysis_name}': 'CAVaLRi_SNPeff_{}'.format(date.today().isoformat()),
                    '{input_vcf}': os.path.join(config['s3_bucket'], config['s3_vcf_input']),
                    '{snpeff_parquet}': config['snpeff_parquet']
                }

            # Iterate through lines of the template and replace keywords with variables
            for line in f:
                line = line.strip()                
                for k, v in map.items():
                    if re.search(k, line):
                        line = line.replace(k, v)

                print(line, file = fo)
                

    # Send the bash script and output definition file
    emr_ip = config['emr_cluster_ip']
    igm_pem = config['igm_pem']
    local_output_def = os.path.join(config['project_root'], config['varhouse_output_definition'])
    emr_bash = os.path.join(config['emr_root'], os.path.basename(config['varhouse_bash']))
    emr_output_def = os.path.join(config['emr_root'], os.path.basename(config['varhouse_output_definition']))
    os.system(f'scp -i {igm_pem} {new_bash} hadoop@{emr_ip}:{emr_bash}')
    os.system(f'scp -i {igm_pem} {local_output_def} hadoop@{emr_ip}:{emr_output_def}')
    
    # Remove current SNPeff data
    os.system('aws s3 rm --recursive --quiet {}'.format(config['snpeff_parquet']))

    # Run SNPeff
    os.system(f'ssh -i {igm_pem} hadoop@{emr_ip} "bash {emr_bash}"')

    print('WAITING ON SNPEFF')
    # Only continue when SNPeff is finished
    while True:
        command = "aws s3 ls {}".format(config['snpeff_parquet'] + '/')
        cmd = shlex.split(command)
        p = subprocess.Popen(cmd,  stdout=subprocess.PIPE)
        out, err = p.communicate()
        if re.search('_SUCCESS', out.decode()):
            break
        time.sleep(60)
        print('STILL WAITING ...')


def start_cluster():
    print('Cluster Started')


def stop_cluster():
    print('Cluster Stopped')
