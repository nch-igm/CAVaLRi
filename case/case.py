# from ..genotype import *
# from ..phenotype import *
import sys
sys.path.append('..')
from genotype import *
from phenotype import *
from .build_case_data import *
from .calculate_moiLR import *
from .calculate_genotypeLR import *
from .calculate_phenotypeLR import *
from .calculate_compositeLR import *
# from .score_phenotypes import *
from .add_rankings import *
from .add_tp import *
from .build_case_summary import *


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
    def __init__(self, cohort, case_id, genotype_path, phenotype_path, genome_build = 'hg38'):
        self.cohort = cohort
        self.case_id = case_id
        self.genome_build = genome_build
        self.genotype = Genotype(self, genotype_path)
        self.phenotype = Phenotype(self, phenotype_path)

    def build_case_data(self):
        self.case_data = build_case_data(self)

    def calculate_moiLR(self):
        self.case_data = calculate_moiLR(self)
    
    def calculate_genotypeLR(self):
        self.case_data = calculate_genotypeLR(self.case_data)

    def calculate_phenotypeLR(self):
        self.case_data = calculate_phenotypeLR(self)

    # def score_phenotypes(self):
    #     self = score_phenotypes(self)

    def calculate_compositeLR(self):
        self.case_data = calculate_compositeLR(self.case_data)
    
    def add_rankings(self):
        self.case_data = add_rankings(self.case_data)

    def add_tp(self):
        self.case_data = add_tp(self)

    def build_case_summary(self):
        self.case_summary = build_case_summary(self.case_data)
        