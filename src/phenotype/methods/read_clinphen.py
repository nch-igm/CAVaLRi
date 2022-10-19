import sys
import pandas as pd

# Import config
sys.path.append('../config')
from config import *

def read_clinphen(case):
    df = pd.read_csv(case.phenotype_path)[['term_name', os.path.join(config['project_root'], config['ic_column'])]]
    df.sort_values(config['ic_column'])
    return df.head(config['hpo_total_upper_bound'])