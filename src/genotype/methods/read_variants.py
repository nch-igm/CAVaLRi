import sys
import os
import re
import pandas as pd
from pkg_resources import parse_version
import vcf
import subprocess
# sys.path.append('../..')
# sys.path.append('..')
# from config import *

def worker(cmd):
    # cmd = f'export PATH={conda_bin}:$PATH && source activate cavalri && {cmd}'
    p = subprocess.Popen(cmd,  stdout=subprocess.PIPE, shell = True, env={'LANGUAGE':'en_US.en', 'LC_ALL':'en_US.UTF-8'})
    p.wait()
    out, err = p.communicate()
    try:
        return out.decode()
    except:
        return err.decode()
    # print(cmd)


def parse_gene(row_info):
    parsed_values = row_info.split(';')
    for pv in parsed_values:
        k,v = pv.split('=')
        if k == 'Gene.refGene':
            return v


def parse_samples(vcf_reader):

    # Initialize result
    res = {}

    # Get the first row and parse samples
    
    for var in vcf_reader:
        for sample in var.samples:

            # Proband
            if re.search('Proband', sample.sample, re.IGNORECASE) or re.search('_P_', sample.sample, re.IGNORECASE):
                res['PROBAND'] = sample.sample

            # Mother
            if re.search('Mother', sample.sample, re.IGNORECASE) or re.search('_M_', sample.sample, re.IGNORECASE):
                res['MOTHER'] = sample.sample

            # Father
            if re.search('Father', sample.sample, re.IGNORECASE) or re.search('_F_', sample.sample, re.IGNORECASE):
                res['FATHER'] = sample.sample

            #TODO Incorporate inputted proband, mother, father

        break

    return(res)



def read_variants(genotype):

    # Normalize the vcf
    normalized_vcf_dir = os.path.join(genotype.case.cohort.root_path, genotype.case.temp_dir, 'normalized_vcfs')
    if not os.path.exists(normalized_vcf_dir):
        os.mkdir(normalized_vcf_dir)
    norm_vcf = os.path.join(normalized_vcf_dir, f'{genotype.case.case_id}.norm.vcf.gz')
    p = worker(f'{genotype.case.bcftools} norm -f {genotype.case.reference_path} -Oz -o {norm_vcf} {genotype.case.genotype.genotype_path}')
    p = worker(f'/igm/apps/htslib/htslib-1.4/bin/tabix {norm_vcf}')

    # Read in vcf
    vcf_reader = vcf.Reader(filename = norm_vcf, compressed=True, encoding='ISO-8859-1')
    # metadata = vcf_reader.metadata

    # Initialize list to capture variants
    var_list = []

    # Get samples in vcf
    samples = parse_samples(vcf_reader)
    
    for var in vcf_reader:
        
        # Variant specific
        start_pos = 3 if re.search('chr', var.CHROM) else 0
        chrom = var.CHROM[start_pos:]
        pos = var.POS
        ref = var.REF
        alt = ','.join([str(i) for i in var.ALT])
        gene = var.INFO['Gene.refGene'][0]
        var_row = [chrom, pos, ref, alt, gene]
        columns = ['CHROM', 'POS', 'REF', 'ALT','GENE']

        # Sample specific
        for sample in samples.values():
            for s in var.samples:
                if s.sample == sample:

                    # Add sample data and add sample name to columns
                    var_row.append({
                        'GT': s.data.GT,
                        'AD': s.data.AD,
                        'AF': round(0 if s.data.AD[1] == 0 else float(s.data.AD[1]/s.data.DP), 3),
                        'DP': s.data.DP
                    })
                    columns.append({v: k for k, v in samples.items()}[sample])

        # Add variant to list, which will be converted back into a data frame
        var_list.append(var_row)
    
    genotype.processed_genotype_path = norm_vcf

    return pd.DataFrame(var_list, columns = columns, dtype=str)
