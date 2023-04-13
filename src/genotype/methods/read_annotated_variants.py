# from contextlib import AsyncExitStack
import re
import sys
import pandas as pd
sys.path.append('.')
from .read_variants import *


def parse_info(row):

    # Parse the info
    info = {}
    info_list = row['INFO'].split(';')
    for i in info_list:
        try:
            info[i.split('=')[0]] = i.split('=')[1]
        except:
            pass

    return info


def parse_gt(row):

    # Parse the genotype
    gt = {}
    format_list = row['FORMAT'].split(':')
    gt_list = row['Proband'].split(':')
    for i, fl in enumerate(format_list):
        gt[fl] = gt_list[i]

    return gt


def read_annotated_variants(input):

    # Read in annotated vcf
    df = read_variants(input)
    df['info'] = df.apply(parse_info, axis = 1)
    df['gt'] = df.apply(parse_gt, axis = 1)

    return df