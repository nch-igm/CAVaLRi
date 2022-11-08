import os
import json
from ...case import Case

def validate_case(inputs):

    path_variables = ['phenotype','vcf','bcftools','conda_bin','reference_path']
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

            validated_input = validate_case(inputs)
            if validated_input != '':
                print(f'Removed case: {case}, {inputs[validated_input]} not a valid path')
            else:
                cs = Case(
                    cohort = cohort, 
                    case_id = case,
                    phenotype_path = inputs['phenotype'],
                    genotype_path = inputs['vcf'],
                    bcftools = inputs['bcftools'],
                    lirical = inputs['lirical'],
                    exomiser = inputs['exomiser'],
                    conda_bin = inputs['conda_bin'],
                    reference_path = inputs['reference_path'],
                    annotations_path = inputs['annotations_path'],
                    biological_sex = inputs['biological_sex'],
                    proband = inputs['proband'],
                    mother = inputs['mother'],
                    father = inputs['father']
                )
                cohort.add_case(cs)
    
    if len(cohort.cases) == 0:
        return False
    return True
