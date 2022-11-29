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
    
    def retrieve_annotated_variants(self):
        self.annotated_variants = retrieve_annotated_variants(self)

    def score_variants(self):
        self.scored_variants = score_variants(self)

    def optimize_scalars(self):
        self.optimal_positions = optimize_scalars(os.path.join(self.output_root, 'summary'))
    
    def get_cases(self):
        if not get_cases(self):
            print('No input files found!')
            sys.exit(1)

    def build_cohort_summary(self):
        self.cohort_summary = build_cohort_summary(self)

    def compile_cohort_variants(self):
        self.cohort_variants = compile_cohort_variants(self)

    def calculate_statistics(self):
        self.stats = calculate_statistics(self)

    def plot_figures(self):
        self.roc_fig, self.roc_data, self.topn_fig, self.topn_data = plot_figures(self)

    def run(self):
        self.get_cases()
        wf = Workflow(self)
        wf.run()
