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
        # with open(os.path.join(self.case.temp_dir, f'{self.case.case_id}.pheno.pickle'), 'wb') as f:
        #     pickle.dump(self.case, f)
        # with open('/igm/home/rsrxs003/rnb/output/BL-280/ba2bf506-ff27-4e28-a0b8-5cc28cd799dc/sample.pheno.pickle', 'wb') as f:
        #     pickle.dump(self.case, f)
        self.scored_phenotypes = score_phenotypes(self.case)

    def calculate_phenotypeLR(self):
        self.phenoLRs = calculate_phenotypeLR(self.case)
