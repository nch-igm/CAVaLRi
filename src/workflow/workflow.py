from .methods import *
import os
import uuid
import subprocess
import pickle
import pandas as pd
import json
import time
import re
import sys

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

    def run(self):
        
        # Set up temporary directory
        # temp_folder = os.path.join(self.cohort.root_path, 'tmp', str(uuid.uuid4()))
        # temp_folder = '/igm/home/rsrxs003/CAVaLRi/tmp/97ba781d-c11e-4711-a7e5-a323167843ef' # IDN_IDV
        # temp_folder = '/igm/home/rsrxs003/CAVaLRi/tmp/4c24ebb0-328c-4c9a-bd51-9970d1d9c289' # deIDN_deIDV
        # temp_folder = '/igm/home/rsrxs003/CAVaLRi/tmp/0bffd2cf-d03a-4368-a03f-23a35125563d'# deIDN_IDV
        temp_folder = '/igm/home/rsrxs003/CAVaLRi/tmp/14264214-7822-4ed7-b828-2a3c3d8770b4' # IDN_deIDV
        
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)

        # Add temp directory to a snakemake yaml file
        snakemake_yaml_path = os.path.join(os.getcwd(), 'src/workflow', 'config.yaml')
        with open(snakemake_yaml_path,'w') as f:
            print(f'tmp_dir: {temp_folder}', file = f)

        # Create case pickle input files
        for case in self.cohort.cases:
            case.temp_dir = temp_folder
            case_plan_path = os.path.join(case.temp_dir, f'{case.case_id}.pickle')
            with open(case_plan_path, 'wb') as f:
                pickle.dump(case, file = f)

        # Run CAVaLRi case pipeline
        # sys.exit(1)
        # res = self.worker(f"/usr/bin/bash {os.path.join(os.getcwd(), 'src/workflow/sm.bash')}")
        
        
        # res = self.worker(f"""
        #     export SGE_ROOT=/igm/apps/sge/sge-8.1.9_install && \
        #     /igm/apps/sge/sge-8.1.9_install/bin/lx-amd64/qsub \
        #         {os.path.join(os.getcwd(), 'src/workflow/sm_qsub.bash')}
        #         """)
        # print(res)
        # job_id = res.split(' ')[2]
        # while True:
        #     res = self.worker(f"""
        #     export SGE_ROOT=/igm/apps/sge/sge-8.1.9_install && \
        #     /igm/apps/sge/sge-8.1.9_install/bin/lx-amd64/qstat
        #     """)
        #     if not re.search(job_id, res):
        #         break
        #     time.sleep(5)


        # Read in populated case data
        for case in self.cohort.cases:
            # if case.case_id in ['DDDP102005', 'DDDP111619', 'DDDP106875', 'DDDP111151', 'DDDP103048', 'DDDP107286', 'DDDP108067', 'DDDP111262', 'DDDP110896', 'DDDP105431', 'DDDP111242', 'DDDP102820', 'DDDP103122', 'DDDP110713', 'DDDP102547', 'DDDP102051', 'DDDP110961', 'DDDP102251', 'DDDP101839', 'DDDP110800', 'DDDP106745', 'DDDP111468', 'DDDP101866', 'DDDP105451', 'DDDP108103', 'DDDP110753', 'DDDP110872', 'DDDP108406', 'DDDP102297', 'DDDP108441', 'DDDP103710', 'DDDP106936', 'DDDP100281', 'DDDP101854', 'DDDP111096', 'DDDP103071', 'DDDP100285', 'DDDP111178', 'DDDP111286', 'DDDP111390', 'DDDP111249', 'DDDP104617', 'DDDP105999', 'DDDP110983', 'DDDP109995', 'DDDP100091', 'DDDP111187', 'DDDP105140', 'DDDP102114', 'DDDP108896', 'DDDP106414', 'DDDP111106', 'DDDP111322', 'DDDP111516', 'DDDP111219', 'DDDP100213', 'DDDP109352', 'DDDP109893', 'DDDP105767', 'DDDP110890', 'DDDP108836', 'DDDP102111', 'DDDP110970', 'DDDP100264', 'DDDP111217', 'DDDP105394', 'DDDP110776', 'DDDP100030', 'DDDP102294', 'DDDP105825', 'DDDP105942', 'DDDP108492', 'DDDP110796', 'DDDP111027', 'DDDP102680', 'DDDP102389', 'DDDP110748', 'DDDP101852', 'DDDP102594', 'DDDP102110', 'DDDP111190', 'DDDP102497', 'DDDP111119', 'DDDP101867', 'DDDP102726', 'DDDP105749', 'DDDP104795', 'DDDP102052', 'DDDP111265', 'DDDP110784', 'DDDP111487', 'DDDP111515', 'DDDP111703', 'DDDP111227', 'DDDP111271', 'DDDP109873', 'DDDP102206', 'DDDP111266', 'DDDP102589', 'DDDP102140', 'DDDP102221', 'DDDP111350', 'DDDP111459', 'DDDP101224', 'DDDP112652', 'DDDP102138', 'DDDP111317', 'DDDP111423', 'DDDP108234', 'DDDP103666', 'DDDP111596', 'DDDP111428', 'DDDP111250', 'DDDP109404', 'DDDP105217', 'DDDP108415', 'DDDP103044', 'DDDP110962', 'DDDP111211', 'DDDP101100', 'DDDP105223', 'DDDP112690', 'DDDP107978', 'DDDP106981', 'DDDP107416', 'DDDP111513', 'DDDP110760', 'DDDP101851', 'DDDP107364', 'DDDP104896', 'DDDP107459', 'DDDP111081', 'DDDP100284', 'DDDP100161', 'DDDP111198', 'DDDP110825']:
            case.temp_folder = temp_folder
            with open(os.path.join(case.temp_folder, f'{case.case_id}.full.pickle'), 'rb') as f:
                case.case_data = pickle.load(f)
            case.case_summary = pd.read_csv(os.path.join(case.temp_folder, f'{case.case_id}.summary.csv'))
        
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
        # self.worker(f'rm -Rf {temp_folder}')
