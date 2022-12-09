from .methods import *
import os
import uuid
import subprocess
import pickle
import pandas as pd
import json
import time
import shlex
import re
import sys

class Workflow:
    """
    This class is meant to house the computational workflow involved with CAVaLRi
    """

    def __init__(self, cohort):
        self.cohort = cohort

    def worker(self, cmd):
        parsed_cmd = shlex.split(cmd)
        cwd = os.path.join(self.cohort.root_path, 'src/workflow')
        p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
        p.wait()
        out, err = p.communicate()
        # return err.decode()
        return out.decode() if out else err.decode()

    def run(self):

        config = self.cohort.config
        
        # Set up temporary directory
        temp_folder = os.path.join(os.path.abspath(os.path.dirname(self.cohort.input_path)), str(uuid.uuid4()))
        
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)

        # Add temp directory to a snakemake yaml file
        snakemake_yaml_path = os.path.join(self.cohort.root_path, 'src/workflow', 'config.yaml')
        with open(snakemake_yaml_path,'w') as f:
            print(f'tmp_dir: {os.path.abspath(temp_folder)}', file = f)

        # Create case pickle input files
        for case in self.cohort.cases:
            case.temp_dir = temp_folder
            case_plan_path = os.path.join(case.temp_dir, f'{case.case_id}.pickle')
            with open(case_plan_path, 'wb') as f:
                pickle.dump(case, file = f)

        # Run CAVaLRi case pipeline
        res = self.worker(f"snakemake --cores {config['cores']} -pk")
        # res = self.worker(f"qsub {os.path.join(self.cohort.root_path, 'src/workflow/sm_qsub.bash')} {config['cores']}")

        # def parse_qstat():
        #     while True:
        #         res = self.worker(f"qstat -r")
        #         if not re.search('cavalri_scheduler', res):
        #             break
        #         time.sleep(10)
                
        
        # time.sleep(10)
        # parse_qstat()

        # Read in populated case data
        for case in self.cohort.cases:
            case.temp_folder = temp_folder
            with open(os.path.join(case.temp_folder, f'{case.case_id}.full.pickle'), 'rb') as f:
                c = pickle.load(f)
            case.case_data = c.case_data
            case.case_summary = c.case_summary
        
        # Run CAVaLRi cohort summary commands
        self.cohort.build_cohort_summary()
        self.cohort.calculate_statistics()
        self.cohort.compile_cohort_variants()
        self.cohort.plot_figures()

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
        self.cohort.roc_fig.savefig(os.path.join(summary_dir, 'roc.png'))
        self.cohort.topn_fig.savefig(os.path.join(summary_dir, 'topn.png'))
        self.cohort.topn_data.to_csv(os.path.join(summary_dir, 'topn_data.csv'), index=False)
            
        # Remove temporary directory
        self.worker(f'rm -Rf {temp_folder}')
        self.worker(f"rm {os.path.join(self.cohort.root_path, 'src/workflow/cavalri/cavalri_scheduler.o*')}")
