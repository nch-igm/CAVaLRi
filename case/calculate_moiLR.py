from dataclasses import dataclass
from itertools import dropwhile
import sys
import os
import pandas as pd
import yaml
import json
import argparse

# Add package locations to sys.path
sys.path.append('.')
sys.path.append('..')


# import local packages
from config import *


def calculate_moiLR(case):

    # Get gender
    gender_df = pd.read_csv(os.path.join(config['project_root'], config['gender_data']))
    if len(gender_df.loc[gender_df['subjectId'] == case.case_id].index) != 0:
        case.gender = gender_df.loc[gender_df['subjectId'] == case.case_id].reset_index(drop=True).loc[0, 'gender']
    else:
        case.gender = 'undefined'

    # Get available parental samples
    first_gene = list(case.case_data['genes'].keys())[0]
    samples = set(['PROBAND', 'MOTHER', 'FATHER']).intersection(set(case.case_data['genes'][first_gene]['variants'][0].keys()))
        
    
    for gene in case.case_data['genes']:
        for d in case.case_data['genes'][gene]['diseases']:

            # Intialize moi LR
            d['moiLR_log10'] = 0

            # Go through each variant
            for i, v in enumerate(case.case_data['genes'][gene]['variants']):

                # Initialize genotype data
                gt_data = {
                    'PROBAND': 'Unavailable',
                    'MOTHER': 'Unavailable',
                    'FATHER': 'Unavailable'
                    }
                
                # Get genotypes
                for gt in samples:
                    gt_data[gt] = v[gt]['GT']

                # Normalize GT by setting all | to /
                for gt in samples:
                    if gt_data[gt].find('|') != -1:
                        g = gt_data[gt].split('|')
                        gt_data[gt] = '/'.join(g)

                # If both parents are present
                if gt_data['MOTHER'] != 'Unavailable' and gt_data['FATHER'] != 'Unavailable':
                    if d['moi'] == 'AD':

                        # De Novo
                        if gt_data['MOTHER'] == '0/0' and gt_data['FATHER'] == '0/0':
                            d['moiLR_log10'] = 1
                            break

                        # Candidate allele present in the parents
                        if gt_data['MOTHER'] != '0/0' or gt_data['FATHER'] != '0/0':
                            d['moiLR_log10'] = -1
            
                    
                    if d['moi'] == 'AR':

                        # De novo
                        # or parents are heterozygous and the proband is homozygous
                        if gt_data['PROBAND'] == '1/1' and gt_data['MOTHER'] in ('0/0','0/1') and gt_data['FATHER'] in ('0/0', '0/1'):
                            d['moiLR_log10'] = 1
                            break

                        # Not homozygous for only variant
                        if gt_data['PROBAND'] != '1/1' and len(case.case_data['genes'][gene]['variants']) == 1:
                            d['moiLR_log10'] = -1
                            break

                        # If parents are homozygous
                        if gt_data['MOTHER'] == '1/1' or gt_data['FATHER'] == '1/1':
                            d['moiLR_log10'] = -1

                        # Check for a compound het
                        compound_het = set()
                        de_novos = 0
                        if gt_data['PROBAND'] == '0/1' and gt_data['MOTHER'] == '0/1':
                            compound_het.add('MOM HIT')
                        if gt_data['PROBAND'] == '0/1' and gt_data['FATHER'] == '0/1':
                            compound_het.add('DAD HIT')
                        if gt_data['PROBAND'] == '0/1' and gt_data['MOTHER'] == '0/0' and gt_data['FATHER'] == '0/0':
                            compound_het.add('DE NOVO')
                            de_novos += 1

                        if i+1 == len(case.case_data['genes'][gene]['variants']): # If variant is the last variant in the gene
                            if len(compound_het) > 1 or de_novos > 1:
                                d['moiLR_log10'] = 1
                            else:
                                d['moiLR_log10'] = -1
                    
                    if d['moi'] == 'XLD':

                        # De novo
                        if gt_data['PROBAND'] in ('0/1', '1/1') and gt_data['MOTHER'] == '0/0' and gt_data['FATHER'] == '0/0':
                            d['moiLR_log10'] = 1
                            break

                        # Variant present in parents
                        if gt_data['MOTHER'] != '0/0' or gt_data['FATHER'] != '0/0':
                            d['moiLR_log10'] = -1

                        
                    if d['moi'] == 'XLR':

                        # Female case
                        if case.gender == 'F':
                            
                            # De novo
                            if gt_data['PROBAND'] == '1/1' and gt_data['MOTHER'] in ('0/0','0/1') and gt_data['FATHER'] in ('0/0', '0/1'):
                                d['moiLR_log10'] = 1
                                break

                            # If parents are heterozygous and the proband is homozygous
                            if gt_data['PROBAND'] == '1/1' and gt_data['MOTHER'] == '0/1' and gt_data['FATHER'] == '0/1':
                                d['moiLR_log10'] = 1
                                break

                            # Not homozygous for only variant
                            if gt_data['PROBAND'] != '1/1' and len(case.case_data['genes'][gene]['variants']) == 1:
                                d['moiLR_log10'] = -1
                                break

                            # If parents are homozygous
                            if gt_data['MOTHER'] == '1/1' or gt_data['FATHER'] == '1/1':
                                d['moiLR_log10'] = -1

                            # Check for a compound het
                            compound_het = set()
                            de_novos = 0
                            if gt_data['PROBAND'] == '0/1' and gt_data['MOTHER'] == '0/1':
                                compound_het.add('MOM HIT')
                            if gt_data['PROBAND'] == '0/1' and gt_data['FATHER'] == '0/1':
                                compound_het.add('DAD HIT')
                            if gt_data['PROBAND'] == '0/1' and gt_data['MOTHER'] == '0/0' and gt_data['FATHER'] == '0/0':
                                compound_het.add('DE NOVO')
                                de_novos += 1

                            if i+1 == len(case.case_data['genes'][gene]['variants']):
                                if len(compound_het) > 1 or de_novos > 1:
                                    d['moiLR_log10'] = 1
                                else:
                                    d['moiLR_log10'] = -1

                            

                        # Male case
                        if case.gender == 'M':

                            # If mother is heterozygous
                            if gt_data['MOTHER'] == '0/1':
                                d['moiLR_log10'] = 1
                                break
                            
                            # If mother is homozygous alt
                            if gt_data['MOTHER'] == '1/1':
                                d['moiLR_log10'] = -1

                

                # If father is missing
                if gt_data['MOTHER'] != 'Unavailable' and gt_data['FATHER'] == 'Unavailable':

                    if d['moi'] == 'AD':

                        # De Novo
                        if gt_data['MOTHER'] == '0/0':
                            d['moiLR_log10'] = .5
                            break

                        # Candidate allele present in the parent
                        if gt_data['MOTHER'] != '0/0':
                            d['moiLR_log10'] = -1
            
                    
                    if d['moi'] == 'AR':

                        # De novo
                        if gt_data['PROBAND'] == '1/1' and gt_data['MOTHER'] == '0/0':
                            d['moiLR_log10'] = 1
                            break

                        # If parent is heterozygous and the proband is homozygous
                        if gt_data['PROBAND'] == '1/1' and gt_data['MOTHER'] == '0/1':
                            d['moiLR_log10'] = .5
                            break

                        # If parent is homozygous
                        if gt_data['MOTHER'] == '1/1':
                            d['moiLR_log10'] = -1

                        # Not homozygous
                        if gt_data['PROBAND'] != '1/1' and len(case.case_data['genes'][gene]['variants']) == 1:
                            d['moiLR_log10'] = -1
                            break

                    
                    if d['moi'] == 'XLD':

                        # De novo
                        if gt_data['PROBAND'] in ('0/1', '1/1') and gt_data['MOTHER'] == '0/0':
                            d['moiLR_log10'] = .5
                            break

                        # Variant present in parent
                        if gt_data['MOTHER'] != '0/0':
                            d['moiLR_log10'] = -1

                        
                    if d['moi'] == 'XLR':

                        # Female case
                        if case.gender == 'F':
                            
                            # If parent is heterozygous and the proband is homozygous
                            if gt_data['PROBAND'] == '1/1' and gt_data['MOTHER'] in ('0/0', '0/1'):
                                d['moiLR_log10'] = .5
                                break

                            # If parent is homozygous
                            if gt_data['MOTHER'] == '1/1':
                                d['moiLR_log10'] = -1


                        # Male case
                        if case.gender == 'M':

                            # If mother is heterozygous
                            if gt_data['MOTHER'] == '0/1':
                                d['moiLR_log10'] = 1
                                break
                            
                            # If mother is homozygous alt
                            if gt_data['MOTHER'] == '1/1':
                                d['moiLR_log10'] = -1
                


                # If mother is missing
                if gt_data['MOTHER'] == 'Unavailable' and gt_data['FATHER'] != 'Unavailable':

                    if d['moi'] == 'AD':

                        # De Novo
                        if gt_data['FATHER'] == '0/0':
                            d['moiLR_log10'] = .5
                            break

                        # Candidate allele present in the parent
                        if gt_data['FATHER'] != '0/0':
                            d['moiLR_log10'] = -1
            
                    
                    if d['moi'] == 'AR':

                        # De novo
                        if gt_data['PROBAND'] == '1/1' and gt_data['FATHER'] == '0/0':
                            d['moiLR_log10'] = 1
                            break

                        # If parent is heterozygous and the proband is homozygous
                        if gt_data['PROBAND'] == '1/1' and gt_data['FATHER'] == '0/1':
                            d['moiLR_log10'] = .5
                            break

                        # Not homozygous with one pathogenic variant
                        if gt_data['PROBAND'] != '1/1' and len(case.case_data['genes'][gene]['variants']) == 1:
                            d['moiLR_log10'] = -1
                            break

                        # If parent is homozygous
                        if gt_data['FATHER'] == '1/1':
                            d['moiLR_log10'] = -1

                    
                    if d['moi'] == 'XLD':

                        # De novo
                        if gt_data['PROBAND'] in ('0/1', '1/1') and gt_data['FATHER'] == '0/0':
                            d['moiLR_log10'] = .5
                            break

                        # Variant present in parent
                        if gt_data['FATHER'] != '0/0':
                            d['moiLR_log10'] = -1

                        
                    if d['moi'] == 'XLR':

                        # Female case
                        if case.gender == 'F':
                            
                            # If parent is heterozygous and the proband is homozygous
                            if gt_data['PROBAND'] == '1/1' and gt_data['FATHER'] in ('0/0', '0/1'):
                                d['moiLR_log10'] = .5
                                break
                            


                        # Male case
                        if case.gender == 'M':

                            # If father is heterozygous
                            if gt_data['FATHER'] == '0/1':
                                d['moiLR_log10'] = -1
                                break
    
    return case.case_data

