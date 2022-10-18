from ..workflow import Workflow
from .methods import *
# # from .hpo_walk import *
import os
# import vcf


class Cohort:
    """
    This class is meant to summarize a cohort consisting of cases inputted by
    the user
    """

    def __init__(self, input_dir, output_dir, config):
        self.cases = set()
        self.input_path = input_dir
        self.output_path = output_dir
        self.root_path = os.getcwd()
        self.config = config

        # Read in big data files
        # 1) HPO
        # 2) HPOA
        # 3) Propogate frequencies up the tree
        

    def add_case(self, case):
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
    
    def get_cases(self):
        if not get_cases(self):
            print('No input files found!')
            sys.exit(1)

    def run(self):
        self.get_cases()
        print([case.case_id for case in self.cases])
        wf = Workflow(self)
        wf.run()
