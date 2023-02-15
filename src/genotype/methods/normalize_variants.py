import os
import re
import shlex
import subprocess


def worker(cmd):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return out.decode() if out else err.decode()


def normalize_variants(genotype):

    config = genotype.case.cohort.config

    # Normalize the vcf
    normalized_vcf_dir = os.path.join(genotype.case.cohort.root_path, genotype.case.temp_dir, 'normalized_vcfs')
    if not os.path.exists(normalized_vcf_dir):
        os.mkdir(normalized_vcf_dir)
    norm_vcf = os.path.join(normalized_vcf_dir, f'{genotype.case.case_id}.norm.vcf.gz')
    p = worker(f"bcftools norm -f {os.path.join(genotype.case.cohort.root_path, config['reference_path'])} -Oz -o {norm_vcf} {genotype.case.genotype.genotype_path}")
    p = worker(f"tabix {norm_vcf}")

    return norm_vcf
