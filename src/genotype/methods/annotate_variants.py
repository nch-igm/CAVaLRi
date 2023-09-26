import os
import re
import shlex
import subprocess
import sys
import json
from numpy import column_stack
import pandas as pd


def worker(cmd):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')


def annotate_variants(genotype):

    config = genotype.case.cohort.config

    # Get genotype path
    genotype_basename = os.path.basename(genotype.genotype_path)
    conda_bin = genotype.case.cohort.conda_bin
    root_path = genotype.case.cohort.root_path

    # Determine output path
    output = os.path.join(genotype.case.temp_dir, f"{genotype.case.case_id}.annotated")
    temp_out = os.path.abspath(f'{output}.hg38_multianno.vcf')
    unzipped_out = os.path.abspath(f'{output}.vcf')

    # Run annovar
    command = f"""
            {os.path.join(conda_bin, 'perl')} \
            {os.path.join(config['annovar_scripts'],'table_annovar.pl')} \
            -vcfinput {genotype.genotype_path} \
            {config['annovar_db']} \
            -buildver {config['genome_build']} \
            --out {output} \
            -remove \
            -protocol refGene,clinvar_20220320,gnomad30_genome \
            -operation g,f,f -nastring . \
            && mv {temp_out} {unzipped_out} \
            && bgzip {unzipped_out} \
            && tabix {unzipped_out}.gz
    """

    p = worker(command)

    return f"{output}.vcf.gz"
    