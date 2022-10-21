import os
import pandas as pd
import re
import vcf

def get_trio_genotype(case, chrom, pos):

    proband, mother, father = 0, 0, 0
    vcf_reader = vcf.Reader(filename = case.genotype.genotype_path)
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
        'proband_gt': res[0],
        'mother_gt': res[1],
        'father_gt': res[2]
    }


def calculate_moiLR(case):

    config = case.cohort.config

    # Intialize result dict
    res = {'subjectId': case.case_id}

    # Read in inheritance data frame
    moi_df = pd.read_csv(os.path.join(case.cohort.root_path, config['moi_db']))
        

    for d in case.case_data['diseases']:

        # Get MOI
        moi = moi_df.loc[moi_df['omimId'] == 'OMIM:' + str(d['omimId'])]
        if moi.empty:
            moi = ''
        else:
            moi = moi.reset_index(drop=True).loc[0, 'moi']
 
        # Intialize gene status dict
        gene = {
            'variants': d['gene_data']['variants'],
            'count': len(d['gene_data']['variants']),
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
            gt_data = get_trio_genotype(case, chrom, int(pos))

            # Normalize GT by setting all | to /
            for gt in ['proband_gt', 'mother_gt', 'father_gt']:
                if chrom == 'chrX':
                    if len(gt_data[gt]) == 1:
                        if gt_data[gt] == 0:
                            gt_data[gt] = '0/1'
                        else:
                            gt_data[gt] = '0/0'
                        gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]
                    else:
                        gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]

                else:
                    gt_data[gt] = gt_data[gt][0] + '/' + gt_data[gt][2]

            # If the patient is homozygous ref
            if gt_data['proband_gt'] == '0/0':
                gene['moiLR_log10'] = -5

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
                    if case.biological_sex == 'F':
                        
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
                    if case.biological_sex == 'M':

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
                    if case.biological_sex == 'F':
                        
                        # If parent is heterozygous and the proband is homozygous
                        if gt_data['proband_gt'] == '1/1' and gt_data['mother_gt'] in ('0/0', '0/1'):
                            gene['moiLR_log10'] = .5
                            break

                        # If parent is homozygous
                        if gt_data['mother_gt'] == '1/1':
                            gene['moiLR_log10'] = -1


                    # Male case
                    if case.biological_sex == 'M':

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
                    if case.biological_sex == 'F':
                        
                        # If parent is heterozygous and the proband is homozygous
                        if gt_data['proband_gt'] == '1/1' and gt_data['father_gt'] in ('0/0', '0/1'):
                            gene['moiLR_log10'] = .5
                            break
                        


                    # Male case
                    if case.biological_sex == 'M':

                        # If father is heterozygous
                        if gt_data['father_gt'] == '0/1':
                            gene['moiLR_log10'] = -1
                            break
        
        # Append disease to result
        res.update({
            d['omimId']: gene['moiLR_log10']
        })
            
    return res
