import re
import pandas as pd
import vcf

def normalize_fam(vcf):
    return vcf


def aggregate_cohort_vcf(cases):
    
    for i, c in enumerate(cases):

        # Create new vcf for the first case
        if i == 0:
            res_vcf = vcf.Reader(filename = c.genotype_path)
        else:
            next_vcf = vcf.Reader(filename = c.genotype_path)
            res_vcf = normalize_fam(next_vcf)


def aggregate_variants(cases):
    
    for case in cases:
        case_variants = case.genotype.variants[['CHROM', 'POS', 'REF', 'ALT']]
        res_df = case_variants if 'res_df' not in locals() else res_df.append(case_variants)

    return res_df.drop_duplicates().reset_index(drop=True)
    