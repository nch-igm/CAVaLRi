import sys
import os
import subprocess

def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')

def test_cohort():
    conda_bin = os.path.join(sys.exec_prefix, 'bin')
    root_path = os.path.join(os.path.dirname(__file__),'..')
    cmd = f"cd {root_path} && {os.path.join(conda_bin, 'python')} cavalri_cohort.py --input_dir {os.path.join(root_path,'example/cohort')} --output_dir {os.path.join(root_path,'example/cohort/output')}"
    worker(cmd)

test_cohort()