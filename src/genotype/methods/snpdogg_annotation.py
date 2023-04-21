import os
import subprocess
import pandas as pd


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')

def get_snpdogg(genotype):

    config = genotype.case.cohort.config
    conda_bin = genotype.case.cohort.conda_bin
    root_path = genotype.case.cohort.root_path
    snpdogg_parquet_path = os.path.join(genotype.case.cohort.root_path, config['snpdogg'])

    # Set up directories
    wrk_dir = os.path.join(genotype.case.temp_dir, 'snpdogg')
    input_folder = os.path.join(wrk_dir, 'input')
    output_folder = os.path.join(wrk_dir, 'output')

    for dir in [wrk_dir, input_folder, output_folder]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # Reduce the vcf to only include splice variants
    vcf_path = os.path.join(genotype.case.temp_dir, f'{genotype.case.case_id}.filtered.vcf.gz')
    snpdogg_input_vcf = os.path.join(input_folder, f"{genotype.case.case_id}.snpdogg.vcf")
    snv_variants = ['missense_variant','nonsynonymous_SNV']
    exp = ' | '.join([f'INFO/ExonicFunc.refGene==\"{x}\"' for x in snv_variants])
    cmd = f"{os.path.join(conda_bin,'bcftools')} filter -i '{exp}' -Ov -o {snpdogg_input_vcf} {vcf_path}"
    p = worker(cmd)

    # Read in the snpdogg parquet
    snpdogg_df = pd.read_parquet(snpdogg_parquet_path, engine = 'pyarrow')

    # Format variants to be annotated
    snv_df = pd.read_csv(snpdogg_input_vcf, sep = '\t', comment = '#').iloc[:,:5]
    snv_df.columns = ['CHROM','POS','ID','REF','ALT']
    snv_df = snv_df.drop(columns = 'ID').astype({'POS':int})
    snv_df['CHROM'] = snv_df['CHROM'].str[3:]

    # Merge dataframes, effectively annotating candidate SNVs
    snpdogg_df = snv_df.merge(snpdogg_df, on = ['CHROM','POS','REF','ALT'], how = 'left').fillna(0) 
    snpdogg_df = snpdogg_df.astype({'POS':str}).rename(columns = {'SCORE':'snpdogg_score'})
    return snpdogg_df
    