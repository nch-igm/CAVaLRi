from ..workflow import Workflow
from .methods import *
import os
import sys


class Cohort:
    """
    This class is meant to summarize a cohort consisting of cases inputted by
    the user
    """

    def __init__(self, input_dir, output_dir, diagnostic_data, config):
        
        self.conda_bin = os.path.join(sys.exec_prefix, 'bin')
        self.cases = set()
        self.input_path = input_dir
        self.output_path = output_dir
        self.diagnostic_data = diagnostic_data
        self.root_path = os.path.dirname(__file__)
        self.root_path = self.root_path[:self.root_path.find('/src')]
        self.config = config

        # Read in big data files
        # 1) HPO
        # 2) HPOA
        # 3) Propogate frequencies up the tree

        # Create MOI df
        hpoa_path = os.path.join(self.root_path, self.config['hpoa'])
        hpoa = pd.read_csv(hpoa_path, sep = '\t', comment = '#').rename(columns = {'database_id':'OMIM'})
        inherit = {
            'HP:0000006':'AD', 
            'HP:0000007':'AR',
            'HP:0001423':'XLD',
            'HP:0001419':'XLR',
            'HP:0001417':'XLD;XLR'
        }

        moi_df = hpoa[hpoa['hpo_id'].isin(inherit.keys())].rename(columns = {'OMIM':'omimId'})[['omimId','hpo_id']]

        def map_inheritence(row, inherit):
            return inherit[row['hpo_id']]

        moi_df['moi'] = moi_df.apply(map_inheritence, inherit=inherit, axis = 1)
        moi_group_df = moi_df[['omimId','moi']].groupby('omimId')['moi'].apply(list).reset_index()
        moi_group_df = moi_group_df[moi_group_df['omimId'].str.contains('OMIM')]

        def join_inherit(row):
            res = ';'.join(list(set(row['moi'])))
            return ';'.join(list(set(res.split(';'))))

        moi_group_df['moi'] = moi_group_df.apply(join_inherit, axis = 1)
        moi_group_df.to_csv(os.path.join(self.root_path,'data','moi.csv'), index = False)
        

    def add_case(self, case):
       self.cases.add(case)
    
    def remove_case(self, case):
        self.cases.remove(case)
    
    # def retrieve_annotated_variants(self):
    #     self.annotated_variants = retrieve_annotated_variants(self)

    # def score_variants(self):
    #     self.scored_variants = score_variants(self)

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
