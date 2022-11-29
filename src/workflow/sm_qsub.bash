#!/bin/bash
#$ -N cavalri_scheduler
#$ -wd /igm/home/rsrxs003/CAVaLRi/src/workflow
#$ -j y

export PATH=$PATH:/igm/home/rsrxs003/miniconda3/condabin/
conda activate cavalri_test

snakemake \
    --cluster "qsub -pe smp 1 -cwd -o qsub_logfiles -e qsub_logfiles" \
    --jobname "sm_{rule}_{wildcards.case}_{jobid}" \
    -j $1 \
    -pk \
    --latency-wait 60