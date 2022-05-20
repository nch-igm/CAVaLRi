import sys
sys.path.append('..')
sys.path.append('.')
import pickle
from cavalri import *

def main(cases):

    # Create cohort
    x = Cohort(
        genome_build='hg38',
        pheno_path = os.path.join(config['project_root'], config['output_root'], config['clinphen_output']),
        vcf_path = os.path.join(config['project_root'], config['vcf_input'])
    )

    # Parse cases
    cases = cases[2:-2].split("', '")

    # Add cases
    for c in cases:
        cs = Case(
            cohort = x, 
            case_id = c,
            phenotype_path = os.path.join(x.pheno_path, c + '.clinphen.tsv'),
            genotype_path = os.path.join(x.vcf_path, c + '.vcf')
        )
        x.add_case(cs)


    return x


if __name__ == "__main__":
    
    # Intialize parser
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--cases', '-c', type=str, help='Cases that should be included in the cohort')
    parser.add_argument('--output', '-o', type=str, help='Path to output cohort pickle')

    args = parser.parse_args()

    x = main(args.cases)

    with open(args.output, 'wb') as f:
        pickle.dump(x, file = f)
