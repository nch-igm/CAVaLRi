from .methods import *
import os
import uuid
import subprocess
import pickle
import pandas as pd
import json

class Workflow:
    """
    This class is meant to house the computational workflow involved with CAVaLRi
    """

    def __init__(self, cohort):
        self.cohort = cohort

    def worker(self, cmd):
        p = subprocess.Popen(cmd,  stdout=subprocess.PIPE, shell = True, env={'LANGUAGE':'en_US.en', 'LC_ALL':'en_US.UTF-8'})
        p.wait()
        out, err = p.communicate()
        try:
            return out.decode()
        except:
            return err.decode()
        # print(cmd)

    def run(self):
        
        # Set up temporary directory
        temp_folder = os.path.join(self.cohort.root_path, str(uuid.uuid4()))
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)

        # Add temp directory to a snakemake yaml file
        snakemake_yaml_path = os.path.join(os.getcwd(), 'src/workflow', 'config.yaml')
        with open(snakemake_yaml_path,'w') as f:
            print(f'tmp_dir: {temp_folder}', file = f)

        # Create case pickle input files
        for case in self.cohort.cases:
            case.temp_dir = temp_folder
            case_plan_path = os.path.join(temp_folder, f'{case.case_id}.pickle')
            with open(case_plan_path, 'wb') as f:
                pickle.dump(case, file = f)

        # Run CAVaLRi case pipeline
        res = self.worker(f"/usr/bin/bash {os.path.join(os.getcwd(), 'src/workflow/sm.bash')}")

        # Read in populated case data
        for case in self.cohort.cases:
            with open(os.path.join(case.temp_dir, f'{case.case_id}.full.pickle'), 'rb') as f:
                case.case_data = pickle.load(f)
            case.case_summary = pd.read_csv(os.path.join(case.temp_dir, f'{case.case_id}.summary.csv'))
        
        # Run CAVaLRi cohort summary commands
        self.cohort.build_cohort_summary()
        self.cohort.calculate_statistics()
        self.cohort.compile_cohort_variants()

        # Create output directories
        cases_dir = os.path.join(self.cohort.output_path, 'cases')
        summary_dir = os.path.join(self.cohort.output_path, 'summary')
        for dir in [cases_dir, summary_dir]:
            if not os.path.exists(dir):
                os.mkdir(dir)

        # Copy output files
        for case in self.cohort.cases:
            with open(os.path.join(cases_dir, f'{case.case_id}.json'), 'w') as f:
                json.dump(case.case_data, f, indent = 4)
            case.case_summary.to_csv(os.path.join(cases_dir, f'{case.case_id}.summary.tsv'), sep = '\t', index = False)
        self.cohort.cohort_summary.to_csv(os.path.join(summary_dir, 'cohort_summary.tsv'), sep = '\t', index = False)
        self.cohort.cohort_variants.to_csv(os.path.join(summary_dir, 'cohort_variants.tsv'), sep = '\t', index = False)
        # self.cohort.figs.roc_fig.save_fig(os.path.join(summary_dir, 'roc.png'))
        # self.cohort.figs.topn_fig.save_fig(os.path.join(summary_dir, 'topn.png'))
            
        # Remove temporary directory
        # self.worker(f'rm -Rf {temp_folder}')
