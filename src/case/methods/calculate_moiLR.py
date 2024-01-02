import os
import pandas as pd
import json
import re
import vcf
import time


# def get_trio_genotype(case, chrom, start, end):

#     proband, mother, father = 0, 0, 0
#     vcf_reader = vcf.Reader(filename = os.path.join(case.cohort.temp_dir, f'{case.case_id}.filtered.vcf.gz'), compressed=True, encoding='ISO-8859-1')
#     for i, s in enumerate(vcf_reader.samples):
#         parent_patterns = {
#             'proband': case.proband,
#             'mother': case.mother,
#             'father': case.father
#             }
#         if re.search(parent_patterns['proband'].lower(), s.lower()):  
#             proband = (i, s)
#         if re.search(parent_patterns['mother'].lower(), s.lower()):
#             mother = (i, s)
#         if re.search(parent_patterns['father'].lower(), s.lower()):
#             father = (i, s)
#     res = []
#     query = vcf_reader.fetch(chrom = chrom, start = start - 1, end = end)
#     for q in query:
#         if (start != end and start == q.POS and end == q.INFO['END']) or start == end:

#             for col in [proband, mother, father]:
#                 if col != 0:
#                     res.append(q.samples[col[0]]['GT'])
#                 else:
#                     res.append('Unavailable')
#             break

#     return {
#         'proband': res[0],
#         'mother': res[1],
#         'father': res[2]
#     }



def calculate_moiLR(case):

    config = case.cohort.config

    # Intialize result dict
    res = {}

    pathogenic_df = case.genotype.pathogenic_variants.copy()
    pathogenic_df['TYPE'] = 'SHORT'
    pathogenic_cnv_df = case.genotype.pathogenic_cnvs.copy()
    pathogenic_cnv_df['TYPE'] = 'CNV'
    pathogenic_cnv_df = pathogenic_cnv_df.rename(columns = {'CHROMOSOME':'CHROM'})
    pathogenic_cnv_df['REF'] = pathogenic_cnv_df['END']
    pathogenic_cnv_df['ALT'] = pathogenic_cnv_df['CHANGE']

    # Determine parent status
    parents = [k for k,v in {'mother':case.mother,'father':case.father}.items() if v != 'Unavailable']
    family = ['proband'] + parents

    # Read in inheritance data frame
    moi_df = case.cohort.moi.copy()
    for gene in case.case_data['genes'].keys():

        # Intialize dict to store count of alternate alleles
        g_short_df = pathogenic_df[pathogenic_df['GENE_ID'] == gene].reset_index(drop=True)
        g_cnv_df = pathogenic_cnv_df[pathogenic_cnv_df['INTERSECTING_GENES'].apply(lambda x: gene in x)].reset_index(drop=True)
        g_cnv_df = g_cnv_df.rename(columns = {'START':'POS','PATHOGENIC':'score'})
        cols = ['CHROM','POS','REF','ALT','proband','score','TYPE']
        cols = cols + parents
        if len(g_short_df.index) > 0 and len(g_cnv_df.index) == 0:
            g_df = g_short_df.copy()
        elif len(g_short_df.index) == 0 and len(g_cnv_df.index) > 0:
            g_df = g_cnv_df.copy()
        elif len(g_short_df.index) > 0 and len(g_cnv_df.index) > 0:
            g_df = pd.concat([g_short_df[cols],g_cnv_df[cols]], ignore_index=True)
        else:
            print('No pathogenic variants')
            sys.exit(1)
        alt_counts = {f"{var['CHROM']}-{var['POS']}-{var['REF']}-{var['ALT']}":{s:0 for s in family} for idx,var in g_df.iterrows()}
        gts = {f"{var['CHROM']}-{var['POS']}-{var['REF']}-{var['ALT']}":{s:json.loads(var[s])['GT'] for s in family} for idx,var in g_df.iterrows()}

        # Go through each variant
        for chrom_pos, s_counts in alt_counts.items():

            chrom, pos1, pos2 = f"chr{chrom_pos.split('-')[0]}", int(chrom_pos.split('-')[1]), chrom_pos.split('-')[2]
            pos2 = pos1 if not pos2.isnumeric() else int(pos2)

            # Get trio genotypes
            # gt_data = get_trio_genotype(case, chrom, pos1, pos2)
            gt_data = { x: gts[chrom_pos][x] for x in family }

            # Normalize GT by setting all | to /, we don't care about phasing for now
            for gt in family:
                
                # Hemizygous state
                if gt_data[gt] in [1,'1']:
                    gt_data[gt] = '1/1'

                # Phased genotype
                if re.search('\|', gt_data[gt]):
                    g = [int(a) if a != '.' else 0 for a in gt_data[gt].split('|')]
                    g.sort()
                    gt_data[gt] = '/'.join([str(a) for a in g])

                # Unphased genotype
                if re.search('/', gt_data[gt]):
                    g = [int(a) if a != '.' else 0 for a in gt_data[gt].split('/')]
                    g.sort()
                    gt_data[gt] = '/'.join([str(a) for a in g])
            
            # Account for genotype edge cases
            for gt in family:
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
                if gt_data[s] == '1/1':
                    s_counts[s] = 2
                elif s == 'proband' and case.biological_sex == 'M' and gt_data[s] in ['0/1','1/0'] and chrom == 'X':
                    s_counts[s] = 2
                elif gt_data[s] == '0/1':
                    s_counts[s] = 1
                
            
        #TODO Add logic to return variants that agree with mode of inheritance if applicable
        proband_alt_count = sum([var['proband'] for var in alt_counts.values()])
        if 'mother' in family:
            mother_alt_count = sum([var['mother'] for var in alt_counts.values()])
        else:
            mother_alt_count = 0
        if 'father' in family:
            father_alt_count = sum([var['father'] for var in alt_counts.values()])
        else:
            father_alt_count = 0

        # Go through each disaese assocaited with the gene
        for d in case.case_data['genes'][gene].keys():
            
            # Get MOI
            moi_ = moi_df.loc[moi_df['omimId'] == 'OMIM:' + str(d)].reset_index(drop=True)
            if moi_.empty:
                moi_ = ''
            else:
                moi_ = moi_.loc[0, 'moi']

            moi = {m: None for m in moi_.split(';')}

            for k in moi_.split(';'):

                # If both parents are present
                if case.mother != 'Unavailable' and case.father != 'Unavailable':

                    if k in ['AD','XLD']:

                        # Parents are unaffected
                        if all([
                            case.mother_affected == 0,
                            case.father_affected == 0,
                            proband_alt_count >= 1,
                            mother_alt_count == 0,
                            father_alt_count == 0,
                        ]):
                            moi[k] = 1

                        # Mother is affected
                        elif all([
                            case.mother_affected == 1,
                            case.father_affected == 0,
                            proband_alt_count >= 1,
                            mother_alt_count >= 1,
                            father_alt_count == 0,
                        ]):
                            moi[k] = 1

                        # Father is affected
                        elif all([
                            case.mother_affected == 0,
                            case.father_affected == 1,
                            proband_alt_count >= 1,
                            mother_alt_count == 0,
                            father_alt_count >= 1
                        ]):
                            moi[k] = 1

                        # Both parents are affected
                        elif all([
                            case.mother_affected == 1,
                            case.father_affected == 1,
                            proband_alt_count >= 1,
                            mother_alt_count >= 1,
                            father_alt_count >= 1
                        ]):
                            moi[k] = 1

                        # Candidate allele present in either of the parents
                        else:
                            moi[k] = -1
            
                    
                    if k == 'AR':

                        # Unaffected parents
                        if all([
                            case.mother_affected == 0,
                            case.father_affected == 0,
                            proband_alt_count >= 2,
                            mother_alt_count <= 1,
                            father_alt_count <= 1
                        ]):
                            moi[k] = 1

                        # Affected mother
                        elif all([
                            case.mother_affected == 1,
                            case.father_affected == 0,
                            proband_alt_count >= 2,
                            mother_alt_count >= 2,
                            father_alt_count <= 1
                        ]):
                            moi[k] = 1

                        # Affected father
                        elif all([
                            case.mother_affected == 0,
                            case.father_affected == 1,
                            proband_alt_count >= 2,
                            mother_alt_count <= 1,
                            father_alt_count >= 2
                        ]):
                            moi[k] = 1

                        # Affected parents
                        elif all([
                            case.mother_affected == 1,
                            case.father_affected == 1,
                            proband_alt_count >= 2,
                            mother_alt_count >= 2,
                            father_alt_count >= 2
                        ]):
                            moi[k] = 1

                            
                        else:
                            moi[k] = -1

                        
                    if k == 'XLR':

                        # Female case
                        if case.biological_sex == 'F':
                            
                            # Unaffected parents
                            if all([
                                case.mother_affected == 0,
                                case.father_affected == 0,
                                proband_alt_count >= 2,
                                mother_alt_count <= 1,
                                father_alt_count == 0
                            ]):
                                moi[k] = 1

                            # Affected mother
                            elif all([
                                case.mother_affected == 1,
                                case.father_affected == 0,
                                proband_alt_count >= 2,
                                mother_alt_count >= 2,
                                father_alt_count == 0
                            ]):
                                moi[k] = 1
                            
                            # Affected father
                            elif all([
                                case.mother_affected == 0,
                                case.father_affected == 1,
                                proband_alt_count >= 2,
                                mother_alt_count <= 1,
                                father_alt_count >= 1
                            ]):
                                moi[k] = 1
                            
                            # Affected parents
                            elif all([
                                case.mother_affected == 1,
                                case.father_affected == 1,
                                proband_alt_count >= 2,
                                mother_alt_count >= 2,
                                father_alt_count >= 1
                            ]):
                                moi[k] = 1
                  
                            else:
                                moi[k] = -1
                            
                        # Male case
                        if case.biological_sex == 'M':

                            # Unaffected parents
                            if all([
                                case.mother_affected == 0,
                                case.father_affected == 0,
                                proband_alt_count >= 1,
                                mother_alt_count <= 1,
                                father_alt_count == 0
                            ]):
                                moi[k] = 1

                            # Affected mother
                            elif all([
                                case.mother_affected == 1,
                                case.father_affected == 0,
                                proband_alt_count >= 1,
                                mother_alt_count >= 2,
                                father_alt_count == 0
                            ]):
                                moi[k] = 1
                            
                            # Affected father
                            elif all([
                                case.mother_affected == 0,
                                case.father_affected == 1,
                                proband_alt_count >= 1,
                                mother_alt_count <= 1,
                                father_alt_count >= 1
                            ]):
                                moi[k] = 1
                            
                            # Affected parents
                            elif all([
                                case.mother_affected == 1,
                                case.father_affected == 1,
                                proband_alt_count >= 1,
                                mother_alt_count >= 2,
                                father_alt_count >= 1
                            ]):
                                moi[k] = 1
                  
                            else:
                                moi[k] = -1                            

                

                # If father is missing
                elif case.mother != 'Unavailable' and case.father == 'Unavailable':

                    if k in ['AD','XLD']:

                        # Mother is unaffected
                        if all([
                            case.mother_affected == 0,
                            proband_alt_count >= 1,
                            mother_alt_count == 0
                        ]):
                            moi[k] = .5

                        # Mother is affected
                        elif all([
                            case.mother_affected == 1,
                            proband_alt_count >= 1,
                            mother_alt_count >= 1
                        ]):
                            moi[k] = .5

                        # Candidate allele present in the parent
                        else:
                            moi[k] = -1
            
                    if k == 'AR':

                        # Mother is unaffected
                        if all([
                            case.mother_affected == 0,
                            proband_alt_count >= 2,
                            mother_alt_count <= 1
                        ]):
                            moi[k] = 0.5
                        
                        # Mother is affected
                        elif all([
                            case.mother_affected == 1,
                            proband_alt_count >= 2,
                            mother_alt_count >= 2
                        ]):
                            moi[k] = .5

                        else:
                            moi[k] = -1

                    if k == 'XLR':

                        # Female case
                        if case.biological_sex == 'F':
                            
                            # Mother is unaffected
                            if all([
                                case.mother_affected == 0,
                                proband_alt_count >= 2,
                                mother_alt_count <= 1
                            ]):
                                moi[k] = 0.5

                            # Mother is affected
                            elif all([
                                case.mother_affected == 1,
                                proband_alt_count >= 2,
                                mother_alt_count >= 2
                            ]):
                                moi[k] = 0.5
                            
                            else:
                                moi[k] = -1
                            
                        # Male case
                        if case.biological_sex == 'M':

                            # Mother is unaffected
                            if all([
                                case.mother_affected == 0,
                                proband_alt_count >= 1,
                                mother_alt_count <= 1
                            ]):
                                moi[k] = 0.5

                            # Mother is affected
                            elif all([
                                case.mother_affected == 1,
                                proband_alt_count >= 1,
                                mother_alt_count >= 2
                            ]):
                                moi[k] = 0.5
                            
                            else:
                                moi[k] = -1
                


                # If mother is missing
                elif case.mother == 'Unavailable' and case.father != 'Unavailable':

                    if k in ['AD','XLD']:

                        # Father is unaffected
                        if all([
                            case.father_affected == 0,
                            proband_alt_count >= 1,
                            father_alt_count == 0
                        ]):
                            moi[k] = .5

                        # Father is affected
                        elif all([
                            case.father_affected == 1,
                            proband_alt_count >= 1,
                            father_alt_count >= 1
                        ]):
                            moi[k] = .5

                        else:
                            moi[k] = -1
            
                    if k == 'AR':

                        # Father is unaffected
                        if all([
                            case.father_affected == 0,
                            proband_alt_count >= 2,
                            father_alt_count <= 1
                        ]):
                            moi[k] = 0.5
                        
                        # Father is affected
                        elif all([
                            case.father_affected == 1,
                            proband_alt_count >= 2,
                            father_alt_count >= 2
                        ]):
                            moi[k] = .5

                        else:
                            moi[k] = -1

                    if k == 'XLR':

                        # Female case
                        if case.biological_sex == 'F':
                            
                            # Father is unaffected
                            if all([
                                case.father_affected == 0,
                                proband_alt_count >= 2,
                                father_alt_count == 0
                            ]):
                                moi[k] = 0.5

                            # Father is affected
                            elif all([
                                case.father_affected == 1,
                                proband_alt_count >= 2,
                                father_alt_count >= 1
                            ]):
                                moi[k] = 0.5
                            
                            else:
                                moi[k] = -1
                            
                        # Male case
                        if case.biological_sex == 'M':

                            # Father is unaffected
                            if all([
                                case.father_affected == 0,
                                proband_alt_count >= 1,
                                father_alt_count == 0
                            ]):
                                moi[k] = 0.5

                            # Father is affected
                            elif all([
                                case.father_affected == 1,
                                proband_alt_count >= 1,
                                father_alt_count >= 1
                            ]):
                                moi[k] = 0.5
                            
                            else:
                                moi[k] = -1
                                
            d_scores = [v for v in moi.values() if v != None]
            res[d] = 0 if len(d_scores) == 0 else max(d_scores)

    case.moiLRs = res
