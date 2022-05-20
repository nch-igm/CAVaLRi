from .read_clinphen import *
from .read_phenotype import *

class Phenotype:    

    # INITIALIZE
    def __init__(self, case, phenotype_path):
        self.phenotype_path = phenotype_path
        self.case = case
    
    def read_phenotypes(self):
        self.phenotypes = read_phenotypes(self.case)
