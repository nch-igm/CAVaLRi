import pandas as pd
import os
import sys
import yaml
import subprocess


def worker(cmd):
    p = subprocess.Popen(cmd,  stdout=subprocess.PIPE, shell = True, env={'LANGUAGE':'en_US.en', 'LC_ALL':'en_US.UTF-8'})
    p.wait()
    out, err = p.communicate()
    try:
        return out.decode()
    except:
        return err.decode()
    # print(cmd)


# def build_yaml(template_path, case_id, vcf, hpo_ids, output_dir, tsv=True):
def build_yaml(case, hpo_ids, tsv=True):

    config = case.cohort.config

    # extract HPO ID#s from the HPO ID column
    with open(os.path.join(case.cohort.root_path, config['lirical_yaml_template'])) as template:
        result = yaml.full_load(template)

    # split template into two pieces
    result1 = dict(result)
    result2 = dict(result)

    # trim each piece
    del result1['hpoIds']
    del result1['negatedHpoIds']
    del result1['prefix']
    del result2['analysis']

    # update "analysis" component
    result1['analysis']['datadir'] = os.path.join(case.cohort.root_path, config['lirical_data_path'])
    result1['analysis']['vcf'] = case.genotype.processed_genotype_path
    result1['analysis']['exomiser'] = config['exomiser_data']

    # update the rest
    result2['hpoIds'] = hpo_ids
    # result2['prefix'] = os.path.join(case.temp_dir, case.case_id)
    result2['prefix'] = case.case_id
    result2['outdir'] = case.temp_dir


    # set tsv parameter
    result1['analysis']['tsv'] = "true" if tsv else "false"

    if tsv:
        yaml_path = os.path.join(case.temp_dir, f'{case.case_id}.tsv.yaml')
        with open(yaml_path,'w') as file:
            yaml.dump(result1, file)
        with open(yaml_path, 'a') as file:
            yaml.dump(result2, file, default_flow_style=None)

    else:
        yaml_path = os.path.join(case.temp_dir, f'{case.case_id}.html.yaml')
        with open(yaml_path,'w') as file:
            yaml.dump(result1, file)
        with open(yaml_path,'a') as file:
            yaml.dump(result2, file, default_flow_style=None)

    return yaml_path


# def run_lirical(case_id, clinphen_df, vcf, hpo_total, output_dir, tsv=True):
def run_lirical(case):

    # Bring in cohort configuration settings
    config = case.cohort.config
    
    # Filter ClinPhen tsv
    filtered_df = case.phenotype.pheno_df

    # Isolate HPO IDs
    hpo_ids = filtered_df['HPO ID'].to_list()

    # Set temp and final filenames
    tsv_temp_filename = os.path.join(case.temp_dir, f'{case.case_id}.tsv')
    tsv_output_filename = os.path.join(case.temp_dir, f'{case.case_id}.lirical.tsv')
    html_temp_filename = os.path.join(case.temp_dir, f'{case.case_id}.html')
    html_output_filename = os.path.join(case.temp_dir, f'{case.case_id}.lirical.html')
    
    try:

        # Build YAML
        tsv_yaml = build_yaml(case = case, hpo_ids = hpo_ids, tsv=True)
        html_yaml = build_yaml(case = case, hpo_ids = hpo_ids, tsv=False)

        for yml in [[tsv_yaml,tsv_temp_filename,tsv_output_filename],
                    [html_yaml,html_temp_filename,html_output_filename]]:

            if not(os.path.exists(yml[2])):

                # Run LIRICAL
                p = worker(f"cd {case.cohort.root_path} && {os.path.join(case.conda_bin,'java')} -Xmx4G -jar {config['lirical_executable']} yaml -y {yml[0]}")

                # Rename results files
                worker(f'cp {yml[1]} {yml[2]}')

        return tsv_output_filename, html_output_filename
    
    except Exception as err:
        print(f'{err}')
        print("Unable to run lirical, file reading error")
        sys.exit(1)




# if __name__ == "__main__":

#     # Intialize parser
#     parser = argparse.ArgumentParser(description='')
#     parser.add_argument('--hpo_list', '-hl', type=str, help='Path to the ordered list of HPO IDs.')
#     parser.add_argument('--output', '-o', type=str, help='Path to a directory where the raw LIRICAL output data will be written.')
#     parser.add_argument('--case_id', '-c', type=str, help='Case ID of the subject')
#     parser.add_argument('--vcf', type=str, help='Path to case vcf')

#     args = parser.parse_args()

#     # Set variables
#     case_id = args.case_id
#     output_dir = args.output
#     hpo_total = config['hpo_total_upper_bound']

#     # Read in ClinPhen results
#     clinphen_df = pd.read_csv(args.hpo_list, sep='\t', header=0)

#     # Run html LIRICAL
#     # run_lirical(case_id=case_id, clinphen_df=clinphen_df, hpo_total=hpo_total, output_dir=output_dir, tsv=False)
#     run_lirical(case_id=case_id, clinphen_df=clinphen_df, vcf=args.vcf, hpo_total=len(clinphen_df.index), output_dir=output_dir, tsv=False)

#     # Run tsv LIRICAL
#     # run_lirical(case_id=case_id, clinphen_df=clinphen_df, hpo_total=hpo_total, output_dir=output_dir, tsv=True)
#     run_lirical(case_id=case_id, clinphen_df=clinphen_df, vcf=args.vcf, hpo_total=len(clinphen_df.index), output_dir=output_dir, tsv=True)