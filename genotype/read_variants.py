import sys
import os
import re
import pandas as pd
import vcf
sys.path.append('../..')
sys.path.append('..')
from config import *

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

        break

    return(res)



def read_variants(vcf_path):

    # Normalize the vcf
    if not os.path.exists(os.path.join(config['project_root'], config['temp_data'], 'normalized_vcfs')):
        os.mkdir(os.path.join(config['project_root'], config['temp_data'], 'normalized_vcfs'))
    norm_vcf = vcf_path[:vcf_path.find('.vcf')] + '.norm.vcf'
    norm_vcf = os.path.join(config['project_root'], config['temp_data'], 'normalized_vcfs', os.path.basename(norm_vcf))
    print('source ~/.bashrc && bcftools norm -f {} {} -o {}'.format(config['human_hg38'], vcf_path, norm_vcf))
    os.system('source ~/.bashrc && bcftools norm -f {} {} -o {}'.format(config['human_hg38'], vcf_path, norm_vcf))

    # Read in vcf
    vcf_reader = vcf.Reader(filename = norm_vcf)
    metadata = vcf_reader.metadata

    # Initialize list to capture variants
    var_list = []

    # Get samples in vcf
    samples = parse_samples(vcf_reader)
    
    for var in vcf_reader:
        
        # Variant specific
        chrom = var.CHROM[3:]
        pos = var.POS
        ref = var.REF
        alt = ','.join([str(i) for i in var.ALT])
        var_row = [chrom, pos, ref, alt]
        columns = ['CHROM', 'POS', 'REF', 'ALT']

        # Sample specific
        for sample in samples.values():
            for s in var.samples:
                if s.sample == sample:

                    # Add sample data and add sample name to columns
                    var_row.append({
                        'GT': s.data.GT,
                        'AD': s.data.AD,
                        'AF': s.data.AF,
                        'DP': s.data.DP
                    })
                    columns.append({v: k for k, v in samples.items()}[sample])

        # Add variant to list, which will be converted back into a data frame
        var_list.append(var_row)

    return pd.DataFrame(var_list, columns = columns, dtype=str)
