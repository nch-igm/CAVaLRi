from .methods import *
import os

class Phenotype:    

    # INITIALIZE
    def __init__(self, case, phenotype_path):
        self.phenotype_path = os.path.abspath(phenotype_path)
        self.case = case

    def read_phenotypes(self):
        self.phenotypes, self.pheno_df = read_phenotypes(self.case)

    def score_phenotypes(self):
        self.scored_phenotypes = score_phenotypes(self)

    def calculate_phenotypeLR(self):
        self.phenoLRs = calculate_phenotypeLR(self.case)
