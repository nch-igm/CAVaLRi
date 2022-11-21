# Check to see if dependency path is provided
dependency_path=${1-dependencies}

# # Download bcftools
# wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 -P $dependency_path
# tar -xjf $dependency_path/bcftools-1.16.tar.bz2 -C $dependency_path
# cd $dependency_path/bcftools-1.16
# ./configure --prefix=$dependency_path
# make
# make install
# cd $dependency_path

# # Download htslib
# wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 -P $dependency_path
# tar -xjf $dependency_path/htslib-1.16.tar.bz2 -C $dependency_path
# cd $dependency_path/htslib-1.16
# ./configure --prefix=$dependency_path
# make
# make install
# cd $dependency_path

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
unzip $lirical_data_path/2209_hg38.zip