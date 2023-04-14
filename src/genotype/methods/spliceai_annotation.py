import os
import glob
import re
import subprocess
import logging
import pandas as pd
import argparse


def run_spliceai_argparser():
    parser = argparse.ArgumentParser(prog=run_spliceai.__name__, description="")
    parser.add_argument('--input', help="S3 path to varhouse output parquet", required=True)
    parser.add_argument('--output', help="S3 path to varhouse output parquet", required=True)
    parser.add_argument('--reference', help="Path to reference", required=True)
    parser.add_argument('--annotation_file', help="Gene annotation file, uses package files when 'grch37' or 'grch38' is passed", required=True)
    return parser

def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')#.replace('\n','')

def run_spliceai(genotype):

    config = genotype.case.cohort.config

    """Run spliceAI predictions pipeline."""
    # parser = run_spliceai_argparser()
    # args = parser.parse_args()
    wrk_dir = os.path.join(genotype.case.temp_dir, 'spliceAI')

    # Set up directories
    input_excel_folder = os.path.join(wrk_dir, 'input')
    vc_folder = os.path.join(wrk_dir, 'predictions')
    resource_folder = os.path.join(wrk_dir, 'resource')
    output_folder = os.path.join(wrk_dir, 'output')

    for dir in [wrk_dir, input_excel_folder, vc_folder, resource_folder, output_folder]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # set up logging
    # LOGGING_FILE = f"{vc_folder}/logging.txt"
    # logging.basicConfig(filename=LOGGING_FILE, level=logging.DEBUG, format="%(asctime)s:%(levelname)s:%(message)s")

    # Download parquet from s3:
    # download_s3_file(args.input_excel, os.path.join(input_excel_folder, 'excel_output.xlsx'))

    # Reduce the vcf to only include splice variants
    vcf_path = os.path.join(genotype.case.temp_dir,f'{genotype.case.case_id}.filtered.vcf.gz')
    spliceai_input_vcf = f"{vcf_path[:vcf_path.find('.vcf')]}.spliceai.vcf"
    cmd = f"{os.path.join(genotype.case.cohort.conda_bin,'bcftools')} filter -i 'INFO/Func.refGene==\"splicing\"' -Ov -o {spliceai_input_vcf} {vcf_path}"
    p = worker(cmd)


    # Run spliceai
    spliceai_output_vcf = f"{vcf_path[:vcf_path.find('.vcf')]}.spliceai_annotated.vcf"
    reference = os.path.join(genotype.case.cohort.root_path,config['reference_path'])
    cmd = f"{os.path.join(genotype.case.cohort.conda_bin,'spliceai')} -I {spliceai_input_vcf} -O {spliceai_output_vcf} -R {reference} -A grch38"
    p = worker(cmd)

    # Read in results
    cols = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    spliceai_df = pd.read_csv(spliceai_output_vcf, sep = '\t', comment = '#')
    spliceai_df.columns = cols + [i for i in range(len(spliceai_df.columns) - len(cols))]
    def parse_info(row):
        info = {x.split('=')[0]:x.split('=')[1] for x in row['INFO'].split(';') if re.search('=',x)}
        try:
            sa = info['SpliceAI'].split('|')
            return max(sa[2:6])
        except:
            return 0

    spliceai_df['spliceai_score'] = spliceai_df.apply(parse_info, axis = 1)
    spliceai_df = spliceai_df[['CHROM','POS','REF','ALT','spliceai_score']].astype({'spliceai_score':float,'POS':str})
    spliceai_df['CHROM'] = spliceai_df['CHROM'].str[3:]
    return spliceai_df

    # Define template vcf
    # template_vcf = '/Users/rsrxs003/projects/LIRICAL/input/EXOME_RERUNS/template.vcf'

    # Convert input from parquet to VCF
    # splicai_vcf_path = excel_to_vcf(args.input_excel, template_vcf, output_folder)

    # download reference fasta file from s3
    # ref_fa = flex_input(args.ref_fasta, reference_folder)
    # ref_dict = '.'.join(args.ref_fasta.split('.')[:-1]) + '.dict'
    # ref_fai = args.ref_fasta + '.fai'

    # flex_input(ref_fai, reference_folder)
    # flex_input(ref_dict, reference_folder)

    # download bam from s3
    # bai_file = '.'.join(args.input_bam.split('.')[:-1]) + '.bai'
    # input_bam = flex_input(args.input_bam, bam_folder)
    # flex_input(bai_file, bam_folder)
    
    # output_vcf = output_vcf
    # spliceai_cmd = f'spliceai -I {splicai_vcf_path} -O {output_folder}/spliceai_scored.vcf -R {args.reference} -A {args.annotation_file}'
    # worker(cmd=spliceai_cmd, warn=True)
    # logging.debug(f"haplotypecaller_cmd: {spliceai_cmd}")

    # list_of_gvcfs = glob.glob(os.path.join(vc_folder, '*g.vcf.gz'))
    # list_of_gvcfs = sorted(list_of_gvcfs)

    # output_gvcf = os.path.join(variants_folder, os.path.basename(args.output_gvcf))
    # output_gvcf_idx = f'{output_gvcf}.tbi'

    # if args.logging:
    #     flex_output(LOGGING_FILE, os.path.dirname(args.logging), os.path.basename(args.logging))
    # flex_output(output_gvcf, os.path.dirname(args.output_gvcf), os.path.basename(args.output_gvcf))
    # flex_output(output_gvcf_idx, os.path.dirname(args.output_gvcf), f'{os.path.basename(args.output_gvcf)}.tbi')

# run_spliceai()