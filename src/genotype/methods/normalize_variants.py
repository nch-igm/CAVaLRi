import os
import re
import shlex
import subprocess


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')


def normalize_variants(genotype):

    config = genotype.case.cohort.config
    conda_bin = genotype.case.cohort.conda_bin

    # Normalize the vcf
    normalized_vcf_dir = os.path.join(genotype.case.temp_dir, 'normalized_vcfs')
    if not os.path.exists(normalized_vcf_dir):
        os.mkdir(normalized_vcf_dir)
    norm_vcf = os.path.join(normalized_vcf_dir, f'{genotype.case.case_id}.norm.vcf.gz')
    p = worker(f"{os.path.join(conda_bin, 'bcftools')} norm -f {config['reference_path']} -Oz -o {norm_vcf} {genotype.genotype_path}")
    p = worker(f"{os.path.join(conda_bin, 'tabix')} {norm_vcf}")

    return norm_vcf
