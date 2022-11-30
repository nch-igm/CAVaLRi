import os
import glob
import logging
import argparse

from s3 import download_s3_file, flex_input, flex_output
from utils import system_exec, excel_to_vcf


def run_spliceai_argparser():
    parser = argparse.ArgumentParser(prog=run_spliceai.__name__, description="")
    parser.add_argument('--input_excel', help="S3 path to varhouse output parquet", required=True)
    parser.add_argument('--reference', help="Path to reference", required=True)
    parser.add_argument('--annotation_file', help="Gene annotation file, uses package files when 'grch37' or 'grch38' is passed", required=True)
    return parser


def run_spliceai():
    """Run spliceAI predictions pipeline."""
    parser = run_spliceai_argparser()
    args = parser.parse_args()
    wrk_dir = os.path.join(os.path.dirname(__file__), 'spliceAI')

    # Set up directories
    input_excel_folder = f"{wrk_dir}/input_excel"
    vc_folder = f"{wrk_dir}/predictions"
    resource_folder = f"{wrk_dir}/resource"
    output_folder = f"{wrk_dir}/resource"
    for dir in [input_excel_folder, vc_folder, resource_folder, output_folder]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # set up logging
    LOGGING_FILE = f"{vc_folder}/logging.txt"
    logging.basicConfig(filename=LOGGING_FILE, level=logging.DEBUG, format="%(asctime)s:%(levelname)s:%(message)s")

    # Download parquet from s3:
    download_s3_file(args.input_excel, os.path.join(input_excel_folder, 'excel_output.xlsx'))

    # Define template vcf
    template_vcf = '/Users/rsrxs003/projects/LIRICAL/input/EXOME_RERUNS/template.vcf'

    # Convert input from parquet to VCF
    splicai_vcf_path = excel_to_vcf(args.input_excel, template_vcf, output_folder)

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
    spliceai_cmd = f'spliceai -I {splicai_vcf_path} -O {output_folder}/spliceai_scored.vcf -R {args.reference} -A {args.annotation_file}'
    system_exec(cmd=spliceai_cmd, warn=True)
    logging.debug(f"haplotypecaller_cmd: {spliceai_cmd}")

    # list_of_gvcfs = glob.glob(os.path.join(vc_folder, '*g.vcf.gz'))
    # list_of_gvcfs = sorted(list_of_gvcfs)

    # output_gvcf = os.path.join(variants_folder, os.path.basename(args.output_gvcf))
    # output_gvcf_idx = f'{output_gvcf}.tbi'

    # if args.logging:
    #     flex_output(LOGGING_FILE, os.path.dirname(args.logging), os.path.basename(args.logging))
    # flex_output(output_gvcf, os.path.dirname(args.output_gvcf), os.path.basename(args.output_gvcf))
    # flex_output(output_gvcf_idx, os.path.dirname(args.output_gvcf), f'{os.path.basename(args.output_gvcf)}.tbi')

run_spliceai()