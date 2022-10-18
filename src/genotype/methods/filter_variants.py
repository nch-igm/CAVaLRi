import sys
import os
import pandas as pd
import re
import vcf

sys.path.append('../../workflow')
from config import *


def on_transcript(row):
    return False if row['SNPEFF_HGVS_C'].find('-') != -1 or row['SNPEFF_HGVS_C'].find('+') != -1 else True

def filter_variants(genotype):

    # Read in annotated variants
    annotated_df = genotype.annotated_variants
    annotated_df = annotated_df.loc[annotated_df['DP'] > 20] # depth filter
    annotated_df = annotated_df.loc[annotated_df['alt_sample_fraction'].astype('float') < 0.05] # igm 50 freq filter
    annotated_df = annotated_df.loc[~(annotated_df['GT'].isin(['0|0', '0/0']))] # proband only filter
    # annotated_df = annotated_df.loc[~((annotated_df['yprob'].astype('float') < 0.6) & (annotated_df['SNPEFF_ANNOTATION'] == 'missense_variant'))] # snpdogg filter
    
    # Only keep rows where gene symbol and transcript agree
    annotated_df['keep'] = annotated_df.apply(on_transcript, axis=1)
    annotated_df = annotated_df.loc[annotated_df['keep']].drop(columns = ['keep']).reset_index(drop=True)

    # Fill NAs
    annotated_df = annotated_df.fillna({
        'yprob': 0.5,
        'gnomad_ex_faf95_popmax': 0
    })

    # Filter on GNOMAD frequency
    annotated_df = annotated_df.loc[annotated_df['gnomad_ex_faf95_popmax'] < config['gnomAD_cutoff']]

    # Filter out non exon variants
    return annotated_df.loc[~annotated_df['SNPEFF_INTERPRETED_EFFECT'].isna()].reset_index(drop=True)
    