import sys
import os
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from genotype.genotype import Genotype
from phenotype.phenotype import Phenotype
from methods import *

class Case:
    """
    Cases correspond to a patient, each of which belong to a 
    parent cohort. Genetic and phenotypic information are captured
    in Genotype and Phenotype classes, respectively. Instances
    of these classes are assigned attributes to a Case class.

    By analyizing the information stored in Genotype and Phenotype 
    classes, disease ranking can be achieved. All data necessary in
    performing this calculation is stored in the Case.case_data attribute.
    """

    # INITIALIZE
    def __init__(self, cohort, case_id, vcf_path, phenotype_path, 
                 biological_sex, proband, mother, father, mother_affected, 
                 father_affected, genome_build='hg38'):
        
        self.cohort = cohort
        self.case_id = case_id
        self.genome_build = genome_build
        self.genotype = Genotype(self, vcf_path)
        self.phenotype = Phenotype(self, phenotype_path)
        self.biological_sex = biological_sex
        self.proband = proband
        self.mother = mother
        self.father = father
        self.mother_affected = int(mother_affected)
        self.father_affected = int(father_affected)

        # Get trio status
        if self.mother != 'Unavailable' and self.father != 'Unavailable':
            self.trio_status = 'TRIO'
        elif (self.mother == 'Unavailable' and self.father != 'Unavailable') | (self.mother != 'Unavailable' and self.father == 'Unavailable'):
            self.trio_status = 'DUO'
        else:
            self.trio_status = 'SINGLETON'

    def validate_inputs(self):
        validate_inputs(self)

    def build_case_data(self):
        build_case_data(self)

    def calculate_moiLR(self):
        calculate_moiLR(self)

    def calculate_compositeLR(self):
        calculate_compositeLR(self)
    
    def add_rankings(self):
        add_rankings(self)

    def build_case_summary(self):
        build_case_summary(self)
        