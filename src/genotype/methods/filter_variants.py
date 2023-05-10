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
        return False if var.samples[proband_pos]['AD'] == None or len(var.samples[proband_pos]['AD']) > 2 else True
    except:
        return True


def igm_common_filter(var, common_ids):
    var_id = f"{var.CHROM}_{var.POS}_{var.REF}_{var.ALT[0]}"
    return True if var_id in common_ids else False


def qual_filter(var, quality_minimum):
    try:
        return True if var.QUAL >= quality_minimum else False
    except:
        return True


def depth_filter(var, depth_minimum, proband):
    try:
        return True if int(var.samples[proband_pos]['DP']) >= depth_minimum else False
    except:
        return True

def synonymous_filter(var):
    return False if var.INFO['ExonicFunc.refGene'][0] == 'synonymous_SNV' else True


def filter_variants(genotype):

    config = genotype.case.cohort.config
    conda_bin = genotype.case.cohort.conda_bin
    root_path = genotype.case.cohort.root_path

    # Set paths
    input_vcf_path = os.path.join(genotype.case.cohort.root_path, genotype.genotype_path)
    filtered_vcf_path = os.path.join(genotype.case.temp_dir, f"{genotype.case.case_id}.filtered.vcf")
    
    # Read in the vcf
    vcf_reader = vcf.Reader(filename = input_vcf_path, compressed=True, encoding='ISO-8859-1')
    vcf_writer = vcf.Writer(open(filtered_vcf_path, 'w'), vcf_reader)

    # Read in IGM variant frequencies
    igm_freq_path = os.path.join(root_path, config['common_variants'])
    if os.path.exists(igm_freq_path):
        igm_common_df = pd.read_csv(igm_freq_path)
        def get_var_id(row):
            return f"chr{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}"
        igm_common_df['var_id'] = igm_common_df.apply(get_var_id, axis = 1)
        igm_common_ids = igm_common_df['var_id'].to_list()
    else:
        igm_common_ids = []

    # Get proband position
    proband_pos = get_proband_pos(vcf_reader.samples, genotype.case.proband)

    # Filter VCF
    for record in vcf_reader:
        if (
            proband_alt_filter(record, proband_pos)
                and
            qual_filter(record, config['quality_minimum'])
                and
            depth_filter(record, config['depth_minimum'], proband_pos)
                and
            multiallelic_filter(record, proband_pos)
                and
                (
                    (
                        exonic_filter(record)
                            and
                        popmax_filter(record, config['gnomAD_popmax'])
                            and
                        synonymous_filter(record)
                            and
                        not igm_common_filter(record, igm_common_ids)
                    )
                        or 
                    clinvar_filter(record)
                )
        ):
            vcf_writer.write_record(record)
    
    # Compress vcf
    vcf_writer.close()

    filtered_vcf_dir = '/igm/home/rsrxs003/rnb/output/BL-283/clinician/filtered_vcfs'
    out_vcf_path = os.path.join(filtered_vcf_dir, f"{genotype.case.case_id}.vcf")
    worker(f"cp {filtered_vcf_path} {filtered_vcf_dir}")

    worker(f"{os.path.join(conda_bin, 'bgzip')} {filtered_vcf_path}")
    worker(f"{os.path.join(conda_bin, 'tabix')} {filtered_vcf_path}.gz")
    return f"{filtered_vcf_path}.gz"
