import sys
import os
import re
import vcf
import subprocess
import pickle


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')


def validate_samples(case):

    # Initialize result
    res = {}
    provided_samples = {k:v for k,v in {'proband':case.proband, 'mother':case.mother, 'father':case.father}.items() if v != 'Unavailable'}

    # Read in vcf
    vcf_reader = vcf.Reader(filename = case.genotype.genotype_path, compressed=True, encoding='ISO-8859-1')

    # Get the first row and parse samples
    
    for var in vcf_reader:
        for sample in var.samples:
            for k,v in provided_samples.items():
                if v == sample.sample:
                    res[k] = sample.sample
        break

    return True if provided_samples == res else False, {k:v for k,v in provided_samples.items() if k not in res.keys()}


def validate_case_paths(case):

    if not os.path.exists(case.phenotype.phenotype_path):
        return 'phenotype'
    if not os.path.exists(case.genotype.genotype_path):
        return 'vcf'
    return ''


def validate_inputs(case):

    # Parse the input file
    output_dir = os.path.dirname(case.cohort.output_path) if not case.cohort.output_path else case.cohort.output_path

    validated_input = validate_case_paths(case)
    samples_pass, samples = validate_samples(case)

    # Report errors
    if validated_input != '':
        print(f"Validation failed: {case.phenotype.phenotype_path}, {validated_input} not a valid path")
        sys.exit(1)

    elif not samples_pass:
        print(f"Validation failed: the following samples were not found in the vcf: {samples}")
        sys.exit(1)

    # Add chr if it chr to CHROM
    conda_bin = case.cohort.conda_bin
    vcf_reader = vcf.Reader(filename = case.genotype.genotype_path, compressed=True, encoding='ISO-8859-1')
    for rec in vcf_reader:
        chrom = rec.CHROM
        break

    if not re.search('chr', chrom):

        # Add chr to a temp file
        o = os.path.join(case.temp_dir, os.path.basename(case.genotype.genotype_path))
        o = o[:o.find('.gz')]
        awk = '{if($0 !~ /^#/) print "chr"$0; else print $0}'
        cmd = f"gunzip -c {case.genotype.genotype_path} | awk '{awk}' > {o}"
        p = worker(cmd)

        # Overwrite the previous
        cmd = f"{os.path.join(conda_bin, 'bgzip')} {o}"
        p = worker(cmd)

    else:

        # Copy the vcf into the temp directory
        o = os.path.join(case.temp_dir, os.path.basename(case.genotype.genotype_path))
        cmd = f"cp {case.genotype.genotype_path} {o}"
        p = worker(cmd)

    return os.path.join(case.temp_dir, os.path.basename(case.genotype.genotype_path))