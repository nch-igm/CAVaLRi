#!/bin/bash
export HOME=/igm/home/rsrxs003
export PATH=$HOME/miniconda3/bin:$PATH
source activate cavalri
cd $HOME/CAVaLRi/src/workflow/
snakemake --cores 1 -pk