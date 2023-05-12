import os
import subprocess
import pandas as pd
import vcf
import re


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')

def get_snv_scores(genotype):

    config = genotype.case.cohort.config
    conda_bin = genotype.case.cohort.conda_bin
    root_path = genotype.case.cohort.root_path
    # snpdogg_parquet_path = os.path.join(genotype.case.cohort.root_path, config['snpdogg'])

    # Set up directories
    wrk_dir = os.path.join(genotype.case.temp_dir, 'snv')
    input_folder = os.path.join(wrk_dir, 'input')
    output_folder = os.path.join(wrk_dir, 'output')

    for dir in [wrk_dir, input_folder, output_folder]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # Reduce the vcf to only include splice variants
    vcf_path = os.path.join(genotype.case.temp_dir, f'{genotype.case.case_id}.filtered.vcf.gz')
    snv_input_vcf = os.path.join(input_folder, f"{genotype.case.case_id}.snv.vcf")
    snv_variants = ['missense_variant','nonsynonymous_SNV']
    exp = ' | '.join([f'INFO/ExonicFunc.refGene==\"{x}\"' for x in snv_variants])
    cmd = f"{os.path.join(conda_bin,'bcftools')} filter -i '{exp}' -Ov -o {snv_input_vcf} {vcf_path}"
    p = worker(cmd)

    # Annotate with dbNSFP
    snv_annotated_vcf = os.path.join(input_folder, f"{genotype.case.case_id}.snv.annotated.vcf")
    cmd = f"""
            {os.path.join(conda_bin, 'perl')} \
            {os.path.join(os.path.join(root_path, config['annovar_scripts']),'table_annovar.pl')} \
            -vcfinput {snv_input_vcf} \
            {os.path.join(root_path, config['annovar_db'])} \
            -buildver {config['genome_build']} \
            --out {snv_annotated_vcf} \
            -remove \
            -protocol dbnsfp42a \
            -operation f -nastring . \
            && mv {snv_annotated_vcf}.{config['genome_build']}_multianno.vcf {snv_annotated_vcf} \
            && {os.path.join(conda_bin, 'bgzip')} {snv_annotated_vcf}
            """
    p = worker(cmd)

    # Format variants to be annotated
    try:
        snv_df = pd.read_csv(snv_input_vcf, sep = '\t', comment = '#', header = None).iloc[:,:5]
        snv_df.columns = ['CHROM','POS','ID','REF','ALT']
        snv_df = snv_df.drop(columns = 'ID').astype({'POS':str})

        # Read in annotated vcf
        vcf_reader = vcf.Reader(filename = f'{snv_annotated_vcf}.gz', compressed=True, encoding='ISO-8859-1')
        var_list = []
        for var in vcf_reader:
            chrom = var.CHROM
            pos = var.POS
            ref = var.REF
            alt = ','.join([str(i) for i in var.ALT])
            try:
                score = var.INFO[config['dbNSFP_score']][0]
                if re.search('|', str(score)):
                    score = str(score).split('|')[0]
            except:
                score = '0'
            var_row = [chrom, pos, ref, alt, score]
            var_list.append(var_row)

        columns = ['CHROM','POS','REF','ALT','SCORE']
        scored_df = pd.DataFrame(var_list, columns = columns, dtype=str)
        scored_df = scored_df.astype({'POS':str})

        # Merge dataframes, effectively annotating candidate SNVs
        scored_snv_df = snv_df.merge(scored_df, on = ['CHROM','POS','REF','ALT'], how = 'left')
        scored_snv_df.loc[scored_snv_df['SCORE'] == 'None','SCORE'] = 0
        scored_snv_df = scored_snv_df.rename(columns = {'SCORE':'snv_score'}).reset_index(drop=True)
        scored_snv_df = scored_snv_df.astype({'snv_score':float,'CHROM':str})
        scored_snv_df['CHROM'] = scored_snv_df['CHROM'].str[3:]
    
    except:
        scored_snv_df = pd.DataFrame()

    return scored_snv_df
    