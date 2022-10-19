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
    def __init__(self, cohort, case_id, genotype_path, phenotype_path, bcftools,
                 lirical, exomiser, conda_bin, reference_path, annotations_path,
                 genome_build='hg38'):

        self.cohort = cohort
        self.case_id = case_id
        self.genome_build = genome_build
        self.bcftools = bcftools
        self.lirical = lirical
        self.exomiser = exomiser
        self.conda_bin = conda_bin
        self.genotype = Genotype(self, genotype_path)
        self.phenotype = Phenotype(self, phenotype_path)
        self.reference_path = reference_path
        self.annotations_path = annotations_path
    
    def run_lirical(self):
        self.lirical_tsv_output, self.lirical_html_output = run_lirical(self)

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
        