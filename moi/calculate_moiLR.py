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
from utils import *


def calculate_moiLR(input):

    # Get json input data
    with open(input, 'r') as lirical_data:
        data = json.load(lirical_data)

    # Set case
    case = data['subjectId']

    # Intialize result dict
    res = {
        'subjectId': case
    }

    # Read in inheritance data frame
    moi_df = pd.read_csv(os.path.join(config['project_root'], config['moi_db']))

    # Get correpsonding vcf
    for entry in os.scandir(config['vcf_ref']):
        vcf_case = entry.name[entry.name.find('_') + 1:entry.name.find('_', entry.name.find('_') + 1)]
        if vcf_case == case and entry.name.find('_P_') == -1 and entry.name.endswith('.vcf'):
            vcf = os.path.join(config['vcf_ref'], entry.name)

    # Get gender
    gender_df = pd.read_csv(os.path.join(config['project_root'], config['gender_data']))
    if len(gender_df.loc[gender_df['subjectId'] == case].index) != 0:
        gender = gender_df.loc[gender_df['subjectId'] == case].reset_index(drop=True).loc[0, 'gender']
    else:
        gender = 'undefined'
    
    for d in data['diseases']:

        # Get MOI
        moi = moi_df.loc[moi_df['omimId'] == 'OMIM:' + str(d['omimId'])]
        if moi.empty:
            moi = ''
        else:
            moi = moi.reset_index(drop=True).loc[0, 'moi']

        # Intialize variant list
        pathogenic_variants = []

        for variant in d['gene_data']['variants']:
            if variant['pathScore'] >= 0.80:
                if moi in ['AD', 'XLD']:
                    if variant['popFreq'] <= config['popFreq_max']:
                        pathogenic_variants.append(variant)
                else:
                    pathogenic_variants.append(variant)

        # Intialize gene status dict
        gene = {
            'variants': pathogenic_variants,
            'count': len(pathogenic_variants),
            'moiLRs': [],
            'moiLR_log10': 0
                        }
        
        # Go through each variant
        for i, pv in enumerate(gene['variants']):

            # Parse out chromosome and position
            chrom_pos = pv['pos']
            chrom = chrom_pos[:chrom_pos.find(':')]
            for j in range(chrom_pos.find(':'), len(chrom_pos)):
                try:
                    x = int(chrom_pos[j])
                    end_pos = j + 1
                except:
                    pass
            pos = chrom_pos[chrom_pos.find(':') + 1: end_pos]

            # Get trio genotypes
            gt_data = get_trio_genotype(vcf, chrom, pos, case)

            # Normalize GT by setting all | to /
            for gt in ['proband_gt', 'mother_gt', 'father_gt']:
                gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]
            
            # If both parents are present
            if gt_data['mother_gt'] != 'U/a' and gt_data['father_gt'] != 'U/a':
                
                if moi == 'AD':

                    # De Novo
                    if gt_data['mother_gt'] == '0/0' and gt_data['father_gt'] == '0/0':
                        gene['moiLR_log10'] = 1
                        break

                    # Candidate allele present in the parents
                    if gt_data['mother_gt'] != '0/0' or gt_data['father_gt'] != '0/0':
                        gene['moiLR_log10'] = -1
        
                
                if moi == 'AR':

                    # De novo
                    if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] in ('0/0','0/1') and gt_data['father_gt'] in ('0/0', '0/1'):
                        gene['moiLR_log10'] = 1
                        break

                    # If parents are heterozygous and the proband is homozygous
                    if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] == '0/1' and gt_data['father_gt'] == '0/1':
                        gene['moiLR_log10'] = 1
                        break

                    # Not homozygous for only variant
                    if gt_data['proband_gt'] != '1/1' and len(gene['variants']) == 1:
                        gene['moiLR_log10'] = -1
                        break

                    # If parents are homozygous
                    if gt_data['mother_gt'] == '1/1' or gt_data['father_gt'] == '1/1':
                        gene['moiLR_log10'] = -1

                    # Check for a compound het
                    if gt_data['proband_gt'] == '0/1' and gt_data['mother_gt'] == '0/1':
                        gene['moiLRs'].append('MOM HIT')
                    if gt_data['proband_gt'] == '0/1' and gt_data['father_gt'] == '0/1':
                        gene['moiLRs'].append('DAD HIT')
                    if gt_data['proband_gt'] == '0/1' and gt_data['mother_gt'] == '0/0' and gt_data['father_gt'] == '0/0':
                        gene['moiLRs'].append('DE NOVO')

                    if i+1 == len(gene['variants']):
                        if len(set(gene['moiLRs'])) > 1 or set(gene['moiLRs']) == set(['DE NOVO']):
                            gene['moiLR_log10'] = 1
                        else:
                            gene['moiLR_log10'] = -1
                
                if moi == 'XLD':

                    # De novo
                    if gt_data['proband_gt'] in ('0/1', '1/1') and gt_data['mother_gt'] == '0/0' and gt_data['father_gt'] == '0/0':
                        gene['moiLR_log10'] = 1
                        break

                    # Variant present in parents
                    if gt_data['mother_gt'] != '0/0' or gt_data['father_gt'] != '0/0':
                        gene['moiLR_log10'] = -1

                    
                if moi == 'XLR':

                    # Female case
                    if gender == 'F':
                        
                        # De novo
                        if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] in ('0/0','0/1') and gt_data['father_gt'] in ('0/0', '0/1'):
                            gene['moiLR_log10'] = 1
                            break

                        # If parents are heterozygous and the proband is homozygous
                        if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] == '0/1' and gt_data['father_gt'] == '0/1':
                            gene['moiLR_log10'] = 1
                            break

                        # Not homozygous for only variant
                        if gt_data['proband_gt'] != '1/1' and len(gene['variants']) == 1:
                            gene['moiLR_log10'] = -1
                            break

                        # If parents are homozygous
                        if gt_data['mother_gt'] == '1/1' or gt_data['father_gt'] == '1/1':
                            gene['moiLR_log10'] = -1

                        # Check for a compound het
                        if gt_data['proband_gt'] == '0/1' and gt_data['mother_gt'] == '0/1':
                            gene['moiLRs'].append('MOM HIT')
                        if gt_data['proband_gt'] == '0/1' and gt_data['father_gt'] == '0/1':
                            gene['moiLRs'].append('DAD HIT')
                        if gt_data['proband_gt'] == '0/1' and gt_data['mother_gt'] == '0/0' and gt_data['father_gt'] == '0/0':
                            gene['moiLRs'].append('DE NOVO')

                        if i+1 == len(gene['variants']):
                            if len(set(gene['moiLRs'])) > 1 or set(gene['moiLRs']) == set(['DE NOVO']):
                                gene['moiLR_log10'] = 1
                            else:
                                gene['moiLR_log10'] = -1

                        

                    # Male case
                    if gender == 'M':

                        # If mother is heterozygous
                        if gt_data['mother_gt'] == '0/1':
                            gene['moiLR_log10'] = 1
                            break
                        
                        # If mother is homozygous alt
                        if gt_data['mother_gt'] == '1/1':
                            gene['moiLR_log10'] = -1

            

            # If father is missing
            if gt_data['mother_gt'] != 'U/a' and gt_data['father_gt'] == 'U/a':

                if moi == 'AD':

                    # De Novo
                    if gt_data['mother_gt'] == '0/0':
                        gene['moiLR_log10'] = .5
                        break

                    # Candidate allele present in the parent
                    if gt_data['mother_gt'] != '0/0':
                        gene['moiLR_log10'] = -1
        
                
                if moi == 'AR':

                    # De novo
                    if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] == '0/0':
                        gene['moiLR_log10'] = 1
                        break

                    # If parent is heterozygous and the proband is homozygous
                    if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] == '0/1':
                        gene['moiLR_log10'] = .5
                        break

                    # If parent is homozygous
                    if gt_data['mother_gt'] == '1/1':
                        gene['moiLR_log10'] = -1

                    # Not homozygous
                    if gt_data['proband_gt'] != '1/1' and len(gene['variants']) == 1:
                        gene['moiLR_log10'] = -1
                        break

                
                if moi == 'XLD':

                    # De novo
                    if gt_data['proband_gt'] in ('0/1', '1/1') and gt_data['mother_gt'] == '0/0':
                        gene['moiLR_log10'] = .5
                        break

                    # Variant present in parent
                    if gt_data['mother_gt'] != '0/0':
                        gene['moiLR_log10'] = -1

                    
                if moi == 'XLR':

                    # Female case
                    if gender == 'F':
                        
                        # If parent is heterozygous and the proband is homozygous
                        if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] in ('0/0', '0/1'):
                            gene['moiLR_log10'] = .5
                            break

                        # If parent is homozygous
                        if gt_data['mother_gt'] == '1/1':
                            gene['moiLR_log10'] = -1


                    # Male case
                    if gender == 'M':

                        # If mother is heterozygous
                        if gt_data['mother_gt'] == '0/1':
                            gene['moiLR_log10'] = 1
                            break
                        
                        # If mother is homozygous alt
                        if gt_data['mother_gt'] == '1/1':
                            gene['moiLR_log10'] = -1
            


            # If mother is missing
            if gt_data['mother_gt'] == 'U/a' and gt_data['father_gt'] != 'U/a':

                if moi == 'AD':

                    # De Novo
                    if gt_data['father_gt'] == '0/0':
                        gene['moiLR_log10'] = .5
                        break

                    # Candidate allele present in the parent
                    if gt_data['father_gt'] != '0/0':
                        gene['moiLR_log10'] = -1
        
                
                if moi == 'AR':

                    # De novo
                    if gt_data['proband_gt'] == '1/1' and gt_data['father_gt'] == '0/0':
                        gene['moiLR_log10'] = 1
                        break

                    # If parent is heterozygous and the proband is homozygous
                    if gt_data['proband_gt'] == '1/1' and gt_data['father_gt'] == '0/1':
                        gene['moiLR_log10'] = .5
                        break

                    # Not homozygous with one pathogenic variant
                    if gt_data['proband_gt'] != '1/1' and len(gene['variants']) == 1:
                        gene['moiLR_log10'] = -1
                        break

                    # If parent is homozygous
                    if gt_data['father_gt'] == '1/1':
                        gene['moiLR_log10'] = -1

                
                if moi == 'XLD':

                    # De novo
                    if gt_data['proband_gt'] in ('0/1', '1/1') and gt_data['father_gt'] == '0/0':
                        gene['moiLR_log10'] = .5
                        break

                    # Variant present in parent
                    if gt_data['father_gt'] != '0/0':
                        gene['moiLR_log10'] = -1

                    
                if moi == 'XLR':

                    # Female case
                    if gender == 'F':
                        
                        # If parent is heterozygous and the proband is homozygous
                        if gt_data['proband_gt'] == '1/1' and gt_data['father_gt'] in ('0/0', '0/1'):
                            gene['moiLR_log10'] = .5
                            break
                        


                    # Male case
                    if gender == 'M':

                        # If father is heterozygous
                        if gt_data['father_gt'] == '0/1':
                            gene['moiLR_log10'] = -1
                            break
        
        # Append disease to result
        res.update({
            d['omimId']: gene['moiLR_log10']
        })
    
    return res



if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input', '-i', type=str, help='Parsed LIRICAL json data')
    parser.add_argument('--output', '-o', type=str, help='Path to MOI appended output json')

    args = parser.parse_args()

    # Get mode of inheritance likelihood ratio (MOI)
    data = calculate_moiLR(args.input)
    with open(args.output, 'w') as json_file:
        json.dump(data, json_file, indent = 4)