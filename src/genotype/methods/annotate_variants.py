import os
import re
import shlex
import subprocess
import sys
import json
from numpy import column_stack
import pandas as pd


def worker(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
    p.wait()
    out, err = p.communicate()
    return err.decode(), out.decode()


def annotate_variants(genotype):

    config = genotype.case.cohort.config

    # Get genotype path
    genotype_basename = os.path.basename(genotype.genotype_path)

    # Determine output path
    output = os.path.join(genotype.case.temp_dir, f"{genotype.case.case_id}.annotated")
    temp_out = os.path.abspath(f'{output}.hg38_multianno.vcf')
    unzipped_out = os.path.abspath(f'{output}.vcf')

    # Run annovar
    command = f"""
        perl {os.path.join(os.path.join(genotype.case.cohort.root_path, config['annovar_scripts']),'table_annovar.pl')} \
            -vcfinput {genotype.genotype_path} \
            {os.path.join(genotype.case.cohort.root_path, config['annovar_db'])} \
            -buildver {config['genome_build']} \
            --out {output} \
            -remove \
            -protocol refGene,clinvar_20220320,gnomad30_genome \
            -operation g,f,f -nastring . \
            && mv {temp_out} {unzipped_out} \
            && bgzip {unzipped_out} \
            && tabix {unzipped_out}.gz
    """
    e,o = worker(command)
    return f"{output}.vcf.gz"
    