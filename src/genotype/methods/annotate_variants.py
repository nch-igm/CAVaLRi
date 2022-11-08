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

# def worker(cmd):
#     p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True, env={'LANGUAGE':'en_US.en', 'LC_ALL':'en_US.UTF-8'})
#     p.wait()
#     out, err = p.communicate()
#     try:
#         return out.decode()
#     except:
#         return err.decode()


# def annovar_annotate_variants(input):

#     # Determine output path
#     output = f"{input[:input.find('.vcf.gz')]}.annotated.vcf"

#     # Run annovar
#     command = f"perl {os.path.join(annovar,'table_annovar.pl')} -vcfinput {input} {os.path.join(annovar,'table_annovar.pl')} humandb/ -buildver hg38 --out {output} -remove -protocol refGene -operation g -nastring ."
#     worker(command)
#     worker(f"mv {output}.hg38_multianno.vcf {output}")
#     return output


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

    # Annotate with ANNOVAR
    # annovar_annotate_variants(genotype.genotype_path)

    # Read in annotated variants
    variant_df = genotype.variants
    ann_df = variant_df.copy()
    # for col in ['gnomad_ex_faf95_popmax','gnomad_wg_faf95_popmax']:
    #     ann_df[col] = ann_df[col].fillna(0)
    # print(ann_df.head(5))

    # ann_df = parse_annotations(genotype.case.annotations_path)

    # # Fill in frequency nulls
    # for col in ['gnomad_ex_faf95_popmax','gnomad_wg_faf95_popmax']:
    #     ann_df[col] = ann_df[col].fillna(0)
    # ann_df = variant_df.merge(ann_df, on = ['CHROM', 'POS', 'REF', 'ALT'])
    # ann_df['DP'] = ann_df.apply(get_depth, axis = 1)
    # ann_df['GT'] = ann_df.apply(get_gt, axis = 1)

    return ann_df
    