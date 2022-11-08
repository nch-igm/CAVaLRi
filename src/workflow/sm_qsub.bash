#!/bin/bash
#$ -N scheduler
#$ -wd /igm/home/rsrxs003/CAVaLRi/src/workflow
#$ -j y

export HOME=/igm/home/rsrxs003
export PATH=$HOME/miniconda3/bin:$PATH
source activate cavalri
snakemake \
    --cluster "qsub -pe smp 1 -cwd -o qsub_logfiles -e qsub_logfiles" \
    --jobname "sm_{rule}_{wildcards.case}_{jobid}" \
    -j 50 \
    -pk \
    --latency-wait 60