import os
import sys
import re

def get_trio_genotype(vcf, chrom, pos, case):

    # Intialize result data
    res = []

    # Read in vcf
    with open(vcf,'r') as f:

        for line in f:
            line=line.rstrip()
            fields=line.split('\t')
            if re.search("^#CHROM",line):
                col_names = line.split('\t')
                samples = col_names[9:]
                proband_col, mother_col, father_col = 0, 0, 0

                # Clinical WES
                # for i, sample in enumerate(samples):
                #     if sample.find('_P_') != -1:
                #         proband_col = i + 9
                #     if sample.find('_M_') != -1:
                #         mother_col = i + 9
                #     if sample.find('_F_') != -1:
                #         father_col = i + 9
                
                # rGS
                for i, sample in enumerate(samples):
                    if sample.find('-Proband-') != -1:
                        proband_col = i + 9
                    if sample.find('-Mother-') != -1:
                        mother_col = i + 9
                    if sample.find('-Father-') != -1:
                        father_col = i + 9

                if proband_col == 0:
                    for i, sample in enumerate(samples):
                        if sample.find('_Proband_') != -1:
                            proband_col = i + 9
                        if sample.find('_Mother_') != -1:
                            mother_col = i + 9
                        if sample.find('_Father_') != -1:
                            father_col = i + 9

            if not re.search("^#", line):

                # Get pertinent fields
                fields=line.split('\t')
                vcf_chrom=fields[0]
                vcf_pos=fields[1]
                
                if vcf_chrom == chrom and pos == vcf_pos:
                    for col in [proband_col, mother_col, father_col]:
                        if col != 0:
                            res.append(fields[col].split(':')[0])
                        else:
                            res.append('Unavailable')
                            
    return {
        'proband_gt': res[0],
        'mother_gt': res[1],
        'father_gt': res[2]
    }
                    
