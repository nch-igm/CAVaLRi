#!/bin/bash
# Check to see if dependency path is provided
dependency_path=${1-dependencies}

# Download reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P $dependency_path
gunzip $dependency_path/hg38.fa.gz
samtools faidx $dependency_path/hg38.fa

# Download annovar databases
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene $dependency_path/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar gnomad30_genome $dependency_path/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar clinvar_20220320 $dependency_path/human_db
dependencies/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar dbnsfp42a $dependency_path/human_db