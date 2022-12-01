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

# Download LIRICAL executable
wget https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/v1.3.4/LIRICAL.jar -P $dependency_path

# Download LIRICAL data
lirical_path=$dependency_path/LIRICAL.jar
lirical_data_path=$dependency_path/lirical_data
java -jar $lirical_path download -d $lirical_data_path

# Download Exomiser data
wget https://data.monarchinitiative.org/exomiser/latest/2209_hg38.zip -P $lirical_data_path
mkdir $lirical_data_path/2209_hg38
unzip $lirical_data_path/2209_hg38.zip -d $lirical_data_path/2209_hg38
rm $lirical_data_path/2209_hg38.zip