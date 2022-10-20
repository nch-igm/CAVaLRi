from .methods import *
import os
import uuid
import subprocess
import pickle

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

        # Create smlinks for input files
        for case in self.cohort.cases:
            case.temp_dir = temp_folder
            case_plan_path = os.path.join(temp_folder, f'{case.case_id}.pickle')
            print(f'Dump path: {os.getcwd()}')
            with open(case_plan_path, 'wb') as f:
                pickle.dump(case, file = f)
                # self.worker(f'ln {e} {case_plan_path}')

        # Copy output files and remove temporary directory
        # print(os.path.join(os.getcwd(), 'src/workflow/sm.bash'))
        res = self.worker(f"/usr/bin/bash {os.path.join(os.getcwd(), 'src/workflow/sm.bash')}")
        # print(res)
        # for case in self.cohort.cases:
        #     self.worker(f"cp {os.path.join(temp_folder, f'{case.case_id}.full.pickle')} {self.cohort.output_path}")
        # self.worker(f'rm -Rf {temp_folder}')
