import os
import vcf
import pickle


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