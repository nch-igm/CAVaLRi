from .aggregate_variants import *
from .aws_commands import *
from .score_variants import *
from.optimize_scalars import *
# from .hpo_walk import *

import sys
sys.path.append('../..')
from config import * 

class Cohort:
    """
    This class is meant to summarize a cohort consisting of cases
    """

    def __init__(self, vcf_path, pheno_path, genome_build = 'hg38'):
        self.cases = set()
        self.vcf_path = vcf_path
        self.pheno_path = pheno_path
        self.genome_build = genome_build
        self.output_root = os.path.join(config['project_root'], config['output_root'])

        # Read in big data files
        # 1) HPO
        # 2) HPOA
        # 3) Propogate frequencies up the tree
        

    def add_case(self, case):
        if case.genome_build != self.genome_build:
            raise 'Incompatible genome build'
        else:
            self.cases.add(case)
    
    def remove_case(self, case):
        self.cases.remove(case)
    
    def aggregate_variants(self):
        self.unique_variants = aggregate_variants(self.cases)
    
    def send_unique_variants(self):
        send_unique_variants(self)
    
    def retrieve_annotated_variants(self):
        self.annotated_variants = retrieve_annotated_variants(self)

    def run_annotation_notebook(self):
        run_annotation_notebook(self)

    def run_snpeff(self):
        run_snpeff(self)

    def score_variants(self):
        self.scored_variants = score_variants(self)

    def optimize_scalars(self):
        self.optimal_positions = optimize_scalars(os.path.join(self.output_root, 'summary'))
        
    def start_cluster(self):
        self.active_cluster = start_cluster()

    def stop_cluster(self):
        stop_cluster(self.active_cluster)
        self.active_cluster = ''
