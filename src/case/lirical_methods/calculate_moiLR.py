import os
import pandas as pd
import re
import vcf

def get_trio_genotype(case, chrom, pos):

    proband, mother, father = 0, 0, 0
    vcf_reader = vcf.Reader(filename = case.genotype.genotype_path, compressed=True, encoding='utf-8')
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
    moiLRs = {d['omimId']:0 for d in case.case_data['diseases']}
    moiLRs['subjectId'] = case.case_id

    # Read in inheritance data frame
    moi_df = pd.read_csv(os.path.join(case.cohort.root_path, config['moi_db']))

    for d in case.case_data['diseases']:

        # Get MOI
        moi = moi_df.loc[moi_df['omimId'] == 'OMIM:' + str(d['omimId'])]
        if moi.empty:
            moi = ''
        else:
            moi = moi.reset_index(drop=True).loc[0, 'moi']
 
        # Intialize dict to store count of alternate alleles
        alt_counts = {var['pos']:{s:0 for s in ['proband', 'mother', 'father']} for var in d['gene_data']['variants']}
        
        # Go through each variant
        for chrom_pos, s_counts in alt_counts.items():

            # Parse out chromosome and position
            chrom = chrom_pos[:chrom_pos.find(':')]
            for j in range(chrom_pos.find(':'), len(chrom_pos)):
                try:
                    x = int(chrom_pos[j])  # Validate that j is an integer
                    end_pos = j + 1
                except:
                    pass
            pos = chrom_pos[chrom_pos.find(':') + 1: end_pos]

            # Get trio genotypes
            gt_data = get_trio_genotype(case, chrom, int(pos))

            # Normalize GT by setting all | to /, we don't care about phasing for now
            for gt in ['proband', 'mother', 'father']:
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

        # If both parents are present
        if case.mother != 'Unavailable' and case.father != 'Unavailable':

            if moi == 'AD':

                # De Novo
                if de_novo_count >= 1:
                    moiLRs[d['omimId']] = 1

                # Candidate allele present in either of the parents
                else:
                    moiLRs[d['omimId']] = -1
    
            
            if moi == 'AR':

                # Check for two hits
                if proband_alt_count >= 2 and mother_alt_count <= 1 and father_alt_count <= 1:
                    moiLRs[d['omimId']] = 1
                else:
                    moiLRs[d['omimId']] = -1
                
            
            if moi == 'XLD':

                # De Novo
                if de_novo_count >= 1:
                    moiLRs[d['omimId']] = 1
                else:
                    moiLRs[d['omimId']] = -1

                
            if moi == 'XLR':

                # Female case
                if case.biological_sex == 'F':
                    
                    # De novo / parents are carriers
                    if proband_alt_count >= 2 and mother_alt_count <= 1 and father_alt_count == 0:
                        moiLRs[d['omimId']] = 1
                    else:
                        moiLRs[d['omimId']] = -1
                     
                # Male case
                if case.biological_sex == 'M':

                    # De Novo or mother is heterozygous
                    if proband_alt_count >= 1 and mother_alt_count <= 1 and father_alt_count == 0:
                        moiLRs[d['omimId']] = 1
                    else:
                        moiLRs[d['omimId']] = -1

        

        # If father is missing
        if case.mother != 'Unavailable' and case.father == 'Unavailable':

            if moi == 'AD':

                # De Novo
                if de_novo_count >= 1:
                    moiLRs[d['omimId']] = .5

                # Candidate allele present in the parent
                else:
                    moiLRs[d['omimId']] = -1
    
            if moi == 'AR':

                # De novo
                if proband_alt_count >= 2 and mother_alt_count <= 1:
                    moiLRs[d['omimId']] = 0.5
                else:
                    moiLRs[d['omimId']] = -1

            if moi == 'XLD':

                # De novo
                if de_novo_count >= 1:
                    moiLRs[d['omimId']] = 0.5
                else:
                    moiLRs[d['omimId']] = -1

            if moi == 'XLR':

                # Female case
                if case.biological_sex == 'F':
                    
                    # If mom has less than 2 hits and the proband has 2 or more hits
                    if proband_alt_count >= 2 and mother_alt_count <= 1:
                        moiLRs[d['omimId']] = 0.5
                    else:
                        moiLRs[d['omimId']] = -1
                    
                # Male case
                if case.biological_sex == 'M':

                    # If mother is heterozygous or de novo
                    if proband_alt_count >= 1 and mother_alt_count <= 1:
                        moiLRs[d['omimId']] = 0.5
                    else:
                        moiLRs[d['omimId']] = -1
        


        # If mother is missing
        if case.mother == 'Unavailable' and case.father != 'Unavailable':

            if moi == 'AD':

                # De Novo
                if de_novo_count >= 1:
                    moiLRs[d['omimId']] = .5

                # Candidate allele present in the parent
                else:
                    moiLRs[d['omimId']] = -1
    
            if moi == 'AR':

                if proband_alt_count >= 2 and father_alt_count <= 1:
                    moiLRs[d['omimId']] = 0.5
                else:
                    moiLRs[d['omimId']] = -1

            if moi == 'XLD':

                # De novo
                if de_novo_count >= 1:
                    moiLRs[d['omimId']] = 0.5
                else:
                    moiLRs[d['omimId']] = -1

            if moi == 'XLR':

                # Female case
                if case.biological_sex == 'F':
                    
                    # If mom has less than 2 hits and the proband has 2 or more hits
                    if proband_alt_count >= 2 and father_alt_count == 0:
                        moiLRs[d['omimId']] = 0.5
                    else:
                        moiLRs[d['omimId']] = -1
                    
                # Male case
                if case.biological_sex == 'M':

                    # If mother is heterozygous or de novo
                    if proband_alt_count >= 1 and father_alt_count == 0:
                        moiLRs[d['omimId']] = 0.5
                    else:
                        moiLRs[d['omimId']] = -1

    return moiLRs
