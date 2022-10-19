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

def annovar_annotate_variants(input):

    # Determine output path
    output = input[:input.find('.vcf')] + '.annotated.vcf'

    # Run annovar
    command = """
		perl {annovar}
		-vcfinput {input}
		{human_db}
		-buildver {build}
		--out {output}
		-remove 
		-protocol refGene
		-operation g
		-nastring . 
    """.format(annovar = config["annovar_script"], input = input, output = output, human_db = config["human_db"], build = config["build"])
    cmd = shlex.split(command)
    p = subprocess.Popen(cmd,  stdout=subprocess.PIPE)
    out, err = p.communicate()
    os.system('mv {} {}'.format(output + '.hg38_multianno.vcf', output))
    return output, out


def parse_annotations(annotations_path):
    try:
        if re.search('.csv', annotations_path):
            return pd.read_csv(annotations_path)
        elif re.search('.xl', annotations_path):
            return pd.read_excel(annotations_path)
        else:
            return pd.read_parquet(annotations_path)

    except Exception as err:
        print(f"Unexpected file type {err=}, {type(err)=}")
        raise


def get_depth(row):
    return json.loads(row['PROBAND'].replace("'", '"'))['DP']

def get_gt(row):
    return json.loads(row['PROBAND'].replace("'", '"'))['GT']

def annotate_variants(genotype):

    # Read in annotated variants
    variant_df = genotype.variants
    ann_df = parse_annotations(genotype.case.annotations_path)

    # Fill in frequency nulls
    for col in ['gnomad_ex_faf95_popmax','gnomad_wg_faf95_popmax']:
        ann_df[col] = ann_df[col].fillna(0)
    ann_df = variant_df.merge(ann_df, on = ['CHROM', 'POS', 'REF', 'ALT'])
    ann_df['DP'] = ann_df.apply(get_depth, axis = 1)
    ann_df['GT'] = ann_df.apply(get_gt, axis = 1)

    return ann_df
    