import os
import re
import shlex
import subprocess
import sys
import json
from numpy import column_stack
import pandas as pd

sys.path.append('../..')
from config import *

def worker(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
    p.wait()
    out, err = p.communicate()
    return err.decode(), out.decode()
    # return out.decode() if out else err.decode()


def annovar_annotate_variants(genotype):

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


# def get_depth(row):
#     return json.loads(row['PROBAND'].replace("'", '"'))['DP']


# def get_gt(row):
#     return json.loads(row['PROBAND'].replace("'", '"'))['GT']


def annotate_variants(genotype):

    # Annotate with ANNOVAR
    new_vcf_path = annovar_annotate_variants(genotype)

    # Read in annotated variants
    # variant_df = genotype.variants
    # ann_df = variant_df.copy()
    # for col in ['gnomad_ex_faf95_popmax','gnomad_wg_faf95_popmax']:
    #     ann_df[col] = ann_df[col].fillna(0)

    # ann_df = parse_annotations(genotype.case.annotations_path)

    # # Fill in frequency nulls
    # for col in ['gnomad_ex_faf95_popmax','gnomad_wg_faf95_popmax']:
    #     ann_df[col] = ann_df[col].fillna(0)
    # ann_df = variant_df.merge(ann_df, on = ['CHROM', 'POS', 'REF', 'ALT'])
    # ann_df['DP'] = ann_df.apply(get_depth, axis = 1)
    # ann_df['GT'] = ann_df.apply(get_gt, axis = 1)

    return new_vcf_path
    