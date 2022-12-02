import os
import subprocess
import argparse
import shlex
import vcf
import pandas as pd
import re

def worker(cmd):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return out.decode() if out else err.decode()


def get_proband_pos(samples, proband):
    for i,s in enumerate(samples):
        if s == proband:
            return i

def proband_alt_filter(var, proband):
    return True if var.samples[proband]['GT'] not in ('0|0','0/0') else False


def popmax_filter(var, thresh):
    pop_max = 0
    for pop in ['AF_afr','AF_ami','AF_amr','AF_asj','AF_eas','AF_fin','AF_nfe','AF_oth','AF_sas']:
        try:
            for p in var.INFO[pop]:
                fp = float(p)
                if pop_max < fp:
                    pop_max = fp
        except:
            pass
    return True if pop_max <= thresh else False


def clinvar_filter(var):
    for clinsig in var.INFO['CLNSIG']:
        if clinsig:
            if re.search('path', clinsig.lower()):
                return True
    return False


def exonic_filter(var):
    for func in var.INFO['Func.refGene']:
        if func in ['exonic','splicing']:
            return True
    return False


def multiallelic_filter(var, proband_pos):
    try:
        return True if var.samples[proband_pos]['AD'] == None or len(var.samples[proband_pos]['AD']) > 2 else False
    except:
        return False


def igm_common_filter(var, common_ids):
    var_id = f"{var.CHROM}_{var.POS}_{var.REF}_{var.ALT[0]}"
    return True if var_id in common_ids else False

def qual_filter(var, quality_minimum):
    try:
        return True if var.QUAL >= quality_minimum else False
    except:
        return True


def filter_variants(genotype):

    config = genotype.case.cohort.config

    # Read in the vcf
    vcf_reader = vcf.Reader(filename = os.path.join(genotype.case.cohort.root_path, genotype.genotype_path), compressed=True, encoding='ISO-8859-1')
    filtered_path = os.path.join(genotype.case.temp_dir, f"{genotype.case.case_id}.filtered.vcf")
    vcf_writer = vcf.Writer(open(filtered_path, 'w'), vcf_reader)


    # Read in IGM variant frequencies
    # igm_common_df = pd.read_csv('/igm/home/rsrxs003/CAVaLRi/data/WES_common_variants.csv')
    # def get_var_id(row):
    #     return f"chr{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}"
    # igm_common_df['var_id'] = igm_common_df.apply(get_var_id, axis = 1)
    # igm_common_ids = igm_common_df['var_id'].to_list()

    # Get proband position
    proband_pos = get_proband_pos(vcf_reader.samples, genotype.case.proband)

    # Filter VCF
    for record in vcf_reader:
        if (
            proband_alt_filter(record, proband_pos)
                and
            qual_filter(record, config['quality_minimum'])
                and
            not multiallelic_filter(record, proband_pos)
                and
                (
                    (
                        exonic_filter(record)
                            and
                        popmax_filter(record, config['gnomAD_popmax'])
                        #     and
                        # not igm_common_filter(record, igm_common_ids)
                    )
                        or 
                    clinvar_filter(record)
                )
        ):
            vcf_writer.write_record(record)
    
    # Compress vcf
    vcf_writer.close()
    worker(f'bgzip {filtered_path}')
    worker(f"tabix {filtered_path}.gz")
    return f"{filtered_path}.gz"
