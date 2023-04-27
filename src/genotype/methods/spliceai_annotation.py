import os
import glob
import re
import subprocess
import logging
import pandas as pd
import argparse


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')


def run_spliceai(genotype):

    config = genotype.case.cohort.config
    conda_bin = genotype.case.cohort.conda_bin
    root_path = genotype.case.cohort.root_path

    """Run spliceAI predictions pipeline."""
    wrk_dir = os.path.join(genotype.case.temp_dir, 'spliceAI')

    # Set up directories
    input_folder = os.path.join(wrk_dir, 'input')
    vc_folder = os.path.join(wrk_dir, 'predictions')
    resource_folder = os.path.join(wrk_dir, 'resource')
    output_folder = os.path.join(wrk_dir, 'output')

    for dir in [wrk_dir, input_folder, vc_folder, resource_folder, output_folder]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # Reduce the vcf to only include splice variants
    vcf_path = os.path.join(genotype.case.temp_dir, f'{genotype.case.case_id}.filtered.vcf.gz')
    spliceai_input_vcf = os.path.join(input_folder, f"{genotype.case.case_id}.spliceai.vcf")
    cmd = f"{os.path.join(conda_bin,'bcftools')} filter -i 'INFO/Func.refGene==\"splicing\"' -Ov -o {spliceai_input_vcf} {vcf_path}"
    p = worker(cmd)

    # Run spliceai
    spliceai_output_vcf = os.path.join(output_folder, f"{genotype.case.case_id}.spliceai_annotated.vcf")
    reference = os.path.join(root_path,config['reference_path'])
    cmd = f"{os.path.join(conda_bin,'spliceai')} -I {spliceai_input_vcf} -O {spliceai_output_vcf} -R {reference} -A grch38"
    p = worker(cmd)

    # Read in results
    cols = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    spliceai_df = pd.read_csv(spliceai_output_vcf, sep = '\t', comment = '#')
    spliceai_df.columns = cols + [i for i in range(len(spliceai_df.columns) - len(cols))]
    
    def parse_info(row):
        try:
            info = {x.split('=')[0]:x.split('=')[1] for x in row['INFO'].split(';') if re.search('=',x)}
            sa = info['SpliceAI'].split('|')
            spliceai_score = max(sa[2:6])
            return 0 if spliceai_score == '.' else spliceai_score
        except:
            return 0

    spliceai_df['spliceai_score'] = spliceai_df.apply(parse_info, axis = 1)
    spliceai_df = spliceai_df[['CHROM','POS','REF','ALT','spliceai_score']].astype({'spliceai_score':float,'POS':str})
    spliceai_df['CHROM'] = spliceai_df['CHROM'].str[3:]
    return spliceai_df
    