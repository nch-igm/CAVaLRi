import os
import shlex
import subprocess
import sys
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


def get_depth(row):
    return row['PROBAND']['DP']

def get_gt(row):
    return row['PROBAND']['GT']

def annotate_variants(genotype):

    # Read in annotated variants
    variant_df = genotype.variants
    annotated_df = variant_df.merge(genotype.case.cohort.annotated_variants, on = ['CHROM', 'POS', 'REF', 'ALT'])
    annotated_df['DP'] = annotated_df.apply(get_depth, axis = 1)
    annotated_df['GT'] = annotated_df.apply(get_gt, axis = 1)

    return annotated_df
    