from .methods import *

class Phenotype:    

    # INITIALIZE
    def __init__(self, case, phenotype_path):
        self.phenotype_path = phenotype_path
        self.case = case

    def read_phenotypes(self):
        self.phenotypes, self.pheno_df = read_phenotypes(self.case)
