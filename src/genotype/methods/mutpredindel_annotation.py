import os
import glob
import re
import subprocess
import logging
import pandas as pd
import argparse
from Bio import SeqIO


def worker(command):
    """
    Runs a bash command using subprocess module
    """
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    return output.decode('utf-8')


def run_mutpredindel(genotype):

    config = genotype.case.cohort.config
    conda_bin = genotype.case.cohort.conda_bin
    root_path = genotype.case.cohort.root_path

    # Set up directories
    wrk_dir = os.path.join(genotype.case.temp_dir, 'mutpredindel')
    input_folder = os.path.join(wrk_dir, 'input')
    vc_folder = os.path.join(wrk_dir, 'predictions')
    resource_folder = os.path.join(wrk_dir, 'resource')
    output_folder = os.path.join(wrk_dir, 'output')

    for dir in [wrk_dir, input_folder, vc_folder, resource_folder, output_folder]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    # Reduce the vcf to only include splice variants
    vcf_path = os.path.join(genotype.case.temp_dir, f'{genotype.case.case_id}.filtered.vcf.gz')
    mutpred_input_vcf = os.path.join(input_folder, f"{genotype.case.case_id}.mutpredindel.vcf")
    inframe_indel_variants = ['nonframeshift_insertion','nonframeshift_deletion',
                        'inframe_insertion','inframe_deletion',
                        'inframe_variant']
    exp = ' | '.join([f'INFO/ExonicFunc.refGene==\"{x}\"' for x in inframe_indel_variants])
    cmd = f"{os.path.join(conda_bin,'bcftools')} filter -i '{exp}' -Ov -o {mutpred_input_vcf} {vcf_path}"
    p = worker(cmd)

    # Read in inframe indel mutations
    inframe_indel_df = pd.read_csv(mutpred_input_vcf, sep = '\t', comment = '#', header = None).iloc[:,[0,1,3,4]]
    inframe_indel_df.columns = ['CHROM','POS','REF','ALT']

    if config['run_mutpredindel']:

        # Create a .avinput file to send to get fasta changes
        annovar_path = os.path.join(root_path, config['annovar_scripts'])
        avinput = os.path.join(input_folder, f"{genotype.case.case_id}.mutpredindel.avinput")
        cmd = f"{os.path.join(annovar_path,'convert2annovar.pl')} -format vcf4 {mutpred_input_vcf} > {avinput}"
        p = worker(cmd)

        # Annotate the .avinput file with exonic functional annotations
        annovar_db = os.path.join(root_path, config['annovar_db'])
        cmd = f"{os.path.join(annovar_path, 'annotate_variation.pl')} -build hg38 {avinput} {annovar_db}"
        p = worker(cmd)

        # Create a fasta input file to send to mutpredindel
        exon_avinput = f"{avinput}.exonic_variant_function"
        refgene_dna = os.path.join(root_path, config['annovar_db'], 'hg38_refGene.txt')
        refgene_mrna = os.path.join(root_path, config['annovar_db'], 'hg38_refGeneMrna.fa')
        mutpred_input_fasta = os.path.join(input_folder, f'{genotype.case.case_id}.input.fasta')
        cmd = f"{os.path.join(annovar_path,'coding_change.pl')} {exon_avinput} {refgene_dna} {refgene_mrna} > {mutpred_input_fasta}"
        p = worker(cmd)

        # Run mutpredindel
        # mutpred_output_vcf = os.path.join(output_folder, f"{genotype.case.case_id}.mutpredindel_annotated.vcf")
        mutpred_script = os.path.join(config['mutpredindel'],'run_MutPredIndel.sh')
        mutpred_mcr = config['muptredindel_MCR']
        output_prefix = os.path.join(output_folder, f'{genotype.case.case_id}.')
        cmd = f"cd {os.path.dirname(mutpred_script)} && {mutpred_script} {mutpred_mcr} {mutpred_input_fasta} {output_prefix}"
        p = worker(cmd)

        # Read in results
        temp_file = os.path.join(output_folder, 'temp_output.txt')
        mutpred_indel_output_path = os.path.join(output_folder,f'{genotype.case.case_id}.output_output.txt')

        with open(temp_file, 'w') as fo:
            with open(mutpred_indel_output_path, 'r') as f:
                for line in f:
                    l = line.split(' ')[:5]
                    info = ''.join(line.split(' ')[5:])
                    score = info.split('|')[1]
                    l.append(score)

        mutpred_score_df = pd.read_csv(temp_file, sep = '\t', header = None)
        mutpred_score_df.columns = ['line','transcript','c_DNA-change','p_change','function','score']

        exon_input_df = pd.read_csv(exonic_variant_input, sep = '\t', header = None).iloc[:,[0,3,4,5,6]]
        exon_input_df.columns = ['line','CHROM','POS','end','type']
        def amend_pos(row):
            if row['type'] != '-':
                return row['start'] - 1
            return row['start']

        mutpred_score_df = exon_input_df.merge(mutpred_score_df, how = 'left').fillna(0)
        mutpred_score_df['start'] = mutpred_score_df.apply(amend_pos, axis = 1)
            
    else:
        
        mutpred_score_df = inframe_indel_df.copy()
        mutpred_score_df['score'] = 0.5
        


    mutpred_score_df = mutpred_score_df[['CHROM','POS','score']]
    mutpred_score_df = inframe_indel_df.merge(mutpred_score_df, on = ['CHROM','POS'])
    mutpred_score_df = mutpred_score_df.rename(columns = {'score':'mutpred_score'})
    mutpred_score_df = mutpred_score_df[['CHROM','POS','REF','ALT','mutpred_score']].astype({'mutpred_score':float,'POS':str})
    mutpred_score_df['CHROM'] = mutpred_score_df['CHROM'].str[3:]
    return mutpred_score_df
    