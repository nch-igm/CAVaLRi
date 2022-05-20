import math

def calculate_compositeLR(case_data):

    for gene in case_data['genes']:
         for d in case_data['genes'][gene]['diseases']:

            # Sum likelihood ratios (allowed since we are in log space)
            d['compositeLR_log10'] = d['moiLR_log10'] + d['phenoLR_log10'] + case_data['genes'][gene]['genotypeLR_log10']

    return case_data
