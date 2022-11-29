import os
import json
from ...case import Case

def validate_case(inputs):

    path_variables = ['phenotype','vcf']
    for pv in path_variables:
        if not os.path.exists(inputs[pv]):
            return pv
    return ''

def get_cases(cohort):

    # Parse cases
    for entry in os.scandir(cohort.input_path):
        if entry.name.endswith('.json'):
            path = os.path.join(cohort.input_path, entry.name)
            case = entry.name[:entry.name.find('.json')]
        
            # Add case
            with open(path,'r') as d:
                inputs = json.load(d)

            # Check for parents
            for parent in ['mother','father']:
                if parent in inputs.keys():
                    if inputs[parent] == '':
                        inputs[parent] = 'Unavailable'
                else:
                    inputs[parent] = 'Unavailable'

            blacklist = ['M17-2917', 'M17-326', 'M19-7871', 'M19-5147', 'M17-2603', 'M19-1819', 'M19-1602', 'M18-5236', 'M18-3212', 'M18-5107', 'M17-3055', 'M20-2011']

            validated_input = validate_case(inputs)
            if validated_input != '':
                print(f'Removed case: {case}, {inputs[validated_input]} not a valid path')

            elif case in blacklist:
                print(f'Removed case: {case}, was on the blacklist')

            else:
                cs = Case(
                    cohort = cohort, 
                    case_id = case,
                    phenotype_path = inputs['phenotype'],
                    genotype_path = inputs['vcf'],
                    biological_sex = inputs['biological_sex'],
                    proband = inputs['proband'],
                    mother = inputs['mother'],
                    father = inputs['father']
                )
                cohort.add_case(cs)
    
    if len(cohort.cases) == 0:
        return False
    return True
