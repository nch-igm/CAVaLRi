import sys
from ..genotype import Genotype
from ..phenotype import Phenotype
from .methods import *


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
    def __init__(self, cohort, case_id, genotype_path, phenotype_path, 
                 biological_sex, proband, mother, father, genome_build='hg38'):
        
        self.cohort = cohort
        self.case_id = case_id
        self.genome_build = genome_build
        self.genotype = Genotype(self, genotype_path)
        self.phenotype = Phenotype(self, phenotype_path)
        self.biological_sex = biological_sex
        self.proband = proband
        self.mother = mother
        self.father = father
    

    def validate_inputs(self):
        self.validation_status = validate_inputs(self)

    def build_case_data(self):
        self.case_data = build_case_data(self)

    def calculate_moiLR(self):
        self.moiLRs = calculate_moiLR(self)

    def calculate_compositeLR(self):
        self.case_data = calculate_compositeLR(self)
    
    def add_rankings(self):
        self.case_data = add_rankings(self)

    def add_tp(self):
        self.case_data = add_tp(self)

    def build_case_summary(self):
        self.case_summary = build_case_summary(self)
        