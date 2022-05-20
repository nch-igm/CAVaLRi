import math

def calculate_genotypeLR(case_data):

    for gene in case_data['genes']:

        # Intialize result
        res = 0.01

        # Iterate through each variant to find the highest score
        for v in case_data['genes'][gene]['variants']:
            if v['max_score'] > res:
                res = v['max_score']

        # Assign genotype LR
        case_data['genes'][gene]['genotypeLR_log10'] = math.log(res, 10)
        
    return case_data
