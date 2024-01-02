import sys
import os
import re
import json
import shlex
import pandas as pd
import vcf
import subprocess


def worker(cmd):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return out.decode() if out else err.decode()


def get_cn_variants(genotype):

    conda_bin = genotype.case.cohort.conda_bin

    genotype.cnv_path = os.path.join(
        genotype.case.cohort.temp_dir, 
        f'{genotype.case.case_id}.cnv.vcf.gz'
        )

    p = worker(f"""
        {os.path.join(conda_bin, 'bcftools')} view -i 'INFO/SVTYPE!=""' -Oz -o {genotype.cnv_path} {genotype.vcf_path}
    """)

    if re.search('Error', p, re.IGNORECASE):
        print(p)
        sys.exit(1)


def get_short_variants(genotype):

    conda_bin = genotype.case.cohort.conda_bin

    genotype.short_variant_path = os.path.join(
        genotype.case.cohort.temp_dir, 
        f'{genotype.case.case_id}.vcf.gz'
        )

    p = worker(f"""
        {os.path.join(conda_bin, 'bcftools')} view -i 'INFO/SVTYPE=""' -Oz -o {genotype.short_variant_path} {genotype.vcf_path}
    """)

    if re.search('Error', p, re.IGNORECASE):
        print(p)
        sys.exit(1)


def split_variants(genotype):

    get_cn_variants(genotype)
    get_short_variants(genotype)