import os
import glob
import logging
import argparse
from utils import system_exec, excel_to_vcf


def run_mutpred_argparser():
    parser = argparse.ArgumentParser(prog=run_mutpred.__name__, description="")
    parser.add_argument('--input', help="Path to a .vcf file containing non-frameshift three nucleotide indels", required=True)
    parser.add_argument('--output', help="Name of output file, in the form of an annotated .vcf", required=True)
    return parser


def run_mutpred():
    """Run spliceAI predictions pipeline."""
    parser = run_mutpred_argparser()
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

    
    # output_vcf = output_vcf
    spliceai_cmd = f'spliceai -I {splicai_vcf_path} -O {output_folder}/mutpred_scored.vcf -R {args.reference} -A {args.annotation_file}'
    system_exec(cmd=spliceai_cmd, warn=True)
    logging.debug(f"haplotypecaller_cmd: {spliceai_cmd}")

   
    # if args.logging:
    #     flex_output(LOGGING_FILE, os.path.dirname(args.logging), os.path.basename(args.logging))
    # flex_output(output_gvcf, os.path.dirname(args.output_gvcf), os.path.basename(args.output_gvcf))
    # flex_output(output_gvcf_idx, os.path.dirname(args.output_gvcf), f'{os.path.basename(args.output_gvcf)}.tbi')

run_mutpred()