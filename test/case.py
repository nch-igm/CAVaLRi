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

def test_case():
    conda_bin = os.path.join(sys.exec_prefix, 'bin')
    root_path = os.path.join(os.path.dirname(__file__),'..')
    cmd = f"{os.path.join(conda_bin, 'python')} {os.path.join(root_path,'cavalri_case.py')} --input {os.path.join(root_path,'example/case/example.json')}"
    print(cmd)

test_case()