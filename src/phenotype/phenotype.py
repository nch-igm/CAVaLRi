from .methods import *
import os
import pickle

class Phenotype:    

    # INITIALIZE
    def __init__(self, case, phenotype_path):
        self.phenotype_path = os.path.join(case.cohort.root_path, phenotype_path)
        self.case = case

    def read_phenotypes(self):
        self.phenotypes, self.pheno_df = read_phenotypes(self.case)

    #TODO score_phenotypes(self.case.genotype.genes, self.phenotypes)
    def score_phenotypes(self):
        self.scored_phenotypes = score_phenotypes(self.case)

    def calculate_phenotypeLR(self):
        self.phenoLRs = calculate_phenotypeLR(self.case)
