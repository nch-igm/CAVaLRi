import os
import pandas as pd
import json
import re
import vcf
import time


def get_trio_genotype(case, chrom, pos):

    proband, mother, father = 0, 0, 0
    vcf_reader = vcf.Reader(filename = os.path.join(case.temp_dir, f'{case.case_id}.filtered.vcf.gz'), compressed=True, encoding='ISO-8859-1')
    for i, s in enumerate(vcf_reader.samples):
        parent_patterns = {
            'proband': case.proband,
            'mother': case.mother,
            'father': case.father
            }
        if re.search(parent_patterns['proband'].lower(), s.lower()):  
            proband = (i, s)
        if re.search(parent_patterns['mother'].lower(), s.lower()):
            mother = (i, s)
        if re.search(parent_patterns['father'].lower(), s.lower()):
            father = (i, s)
    res = []
    query = vcf_reader.fetch(chrom = chrom, start = pos - 1, end = pos)
    for q in query:
        for col in [proband, mother, father]:
            if col != 0:
                res.append(q.samples[col[0]]['GT'])
            else:
                res.append('Unavailable')
        break

    return {
        'proband': res[0],
        'mother': res[1],
        'father': res[2]
    }


def calculate_moiLR(case):

    config = case.cohort.config

    # Intialize result dict
    res = {}

    # Read in inheritance data frame
    moi_df = pd.read_csv(os.path.join(case.cohort.root_path, config['moi_db']))
    for gene in case.case_data['genes'].keys():

        # Intialize dict to store count of alternate alleles
        g_df = case.genotype.variants[case.genotype.variants['GENE'] == gene]
        alt_counts = {f"{var['CHROM']}-{var['POS']}-{var['REF']}-{var['ALT']}":{s:0 for s in ['proband', 'mother', 'father']} for idx,var in g_df.iterrows()}

        # Go through each variant
        for chrom_pos, s_counts in alt_counts.items():

            chrom, pos = f"chr{chrom_pos.split('-')[0]}", chrom_pos.split('-')[1]

            # Get trio genotypes
            gt_data = get_trio_genotype(case, chrom, int(pos))

            # Normalize GT by setting all | to /, we don't care about phasing for now
            for gt in ['proband','mother','father']:
                if gt_data[gt].find('|') != -1:
                    g = gt_data[gt].split('|')
                    gt_data[gt] = '/'.join(g)
                if gt_data[gt].find('/') != -1:
                    g = [int(a) if a != '.' else 0 for a in gt_data[gt].split('/')]
                    g.sort()
                    gt_data[gt] = '/'.join([str(a) for a in g])
            
            # Account for genotype edge cases
            for gt in ['proband','mother','father']:
                if chrom == 'chrX' or chrom == 'X':
                    if len(gt_data[gt]) == 1:
                        if gt_data[gt] == 0:
                            gt_data[gt] = '0/1'
                        else:
                            gt_data[gt] = '0/0'
                        gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]
                    else:
                        gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]

                else:
                    if gt_data[gt] in ['.','1']:
                        gt_data[gt] = '0/1'
                    else:
                        gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]

            # Count alternate alleles
            for s in s_counts.keys():
                if gt_data[s] == '0/1':
                    s_counts[s] = 1
                if gt_data[s] == '1/1':
                    s_counts[s] = 2
            
        #TODO Add logic to return variants that agree with mode of inheritance if applicable
        proband_alt_count = sum([var['proband'] for var in alt_counts.values()])
        mother_alt_count = sum([var['mother'] for var in alt_counts.values()])
        father_alt_count = sum([var['father'] for var in alt_counts.values()])
        de_novo_count = sum([var['proband'] for var in alt_counts.values() if var['mother'] == 0 and var['father'] == 0])

        # Go through each disaese assocaited with the gene
        for d in case.case_data['genes'][gene].keys():
            
            # Get MOI
            moi = moi_df.loc[moi_df['omimId'] == 'OMIM:' + str(d)].reset_index(drop=True)
            if moi.empty:
                moi = ''
            else:
                moi = moi.loc[0, 'moi']

            # If both parents are present
            if case.mother != 'Unavailable' and case.father != 'Unavailable':

                if moi == 'AD':

                    # De Novo
                    if de_novo_count >= 1:
                        res[d] = 1

                    # Candidate allele present in either of the parents
                    else:
                        res[d] = -1
        
                
                if moi == 'AR':

                    # Check for two hits
                    if proband_alt_count >= 2 and mother_alt_count <= 1 and father_alt_count <= 1:
                        res[d] = 1
                    else:
                        res[d] = -1
                    
                
                if moi == 'XLD':

                    # De Novo
                    if de_novo_count >= 1:
                        res[d] = 1
                    else:
                        res[d] = -1

                    
                if moi == 'XLR':

                    # Female case
                    if case.biological_sex == 'F':
                        
                        # De novo / parents are carriers
                        if proband_alt_count >= 2 and mother_alt_count <= 1 and father_alt_count == 0:
                            res[d] = 1
                        else:
                            res[d] = -1
                        
                    # Male case
                    if case.biological_sex == 'M':

                        # De Novo or mother is heterozygous
                        if proband_alt_count >= 1 and mother_alt_count <= 1 and father_alt_count == 0:
                            res[d] = 1
                        else:
                            res[d] = -1

            

            # If father is missing
            if case.mother != 'Unavailable' and case.father == 'Unavailable':

                if moi == 'AD':

                    # De Novo
                    if de_novo_count >= 1:
                        res[d] = .5

                    # Candidate allele present in the parent
                    else:
                        res[d] = -1
        
                if moi == 'AR':

                    # De novo
                    if proband_alt_count >= 2 and mother_alt_count <= 1:
                        res[d] = 0.5
                    else:
                        res[d] = -1

                if moi == 'XLD':

                    # De novo
                    if de_novo_count >= 1:
                        res[d] = 0.5
                    else:
                        res[d] = -1

                if moi == 'XLR':

                    # Female case
                    if case.biological_sex == 'F':
                        
                        # If mom has less than 2 hits and the proband has 2 or more hits
                        if proband_alt_count >= 2 and mother_alt_count <= 1:
                            res[d] = 0.5
                        else:
                            res[d] = -1
                        
                    # Male case
                    if case.biological_sex == 'M':

                        # If mother is heterozygous or de novo
                        if proband_alt_count >= 1 and mother_alt_count <= 1:
                            res[d] = 0.5
                        else:
                            res[d] = -1
            


            # If mother is missing
            if case.mother == 'Unavailable' and case.father != 'Unavailable':

                if moi == 'AD':

                    # De Novo
                    if de_novo_count >= 1:
                        res[d] = .5

                    # Candidate allele present in the parent
                    else:
                        res[d] = -1
        
                if moi == 'AR':

                    if proband_alt_count >= 2 and father_alt_count <= 1:
                        res[d] = 0.5
                    else:
                        res[d] = -1

                if moi == 'XLD':

                    # De novo
                    if de_novo_count >= 1:
                        res[d] = 0.5
                    else:
                        res[d] = -1

                if moi == 'XLR':

                    # Female case
                    if case.biological_sex == 'F':
                        
                        # If mom has less than 2 hits and the proband has 2 or more hits
                        if proband_alt_count >= 2 and father_alt_count == 0:
                            res[d] = 0.5
                        else:
                            res[d] = -1
                        
                    # Male case
                    if case.biological_sex == 'M':

                        # If mother is heterozygous or de novo
                        if proband_alt_count >= 1 and father_alt_count == 0:
                            res[d] = 0.5
                        else:
                            res[d] = -1
            res[d] = 0

    return res
