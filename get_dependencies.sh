#!/bin/bash

# Initialize optional variables with default values
output_dir=""
humandb=""
config_file="$(dirname "$(readlink -f "$0")")/config/config.yaml"
data_dir="$(dirname "$(readlink -f "$0")")/data"


while getopts "a:o:h:" opt; do
  case $opt in
    a)
      annovar="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    h)
      humandb="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if the required option (-a) is provided
if [ -z "$annovar" ]; then
  echo "Error: The -a option (annovar) is required."
  exit 1
fi

# Set default values for output_dir and humandb if not provided
if [ -z "$output_dir" ]; then
  # Get the directory of the script using $0 and then join with 'dependencies'
  script_dir="$(dirname "$(readlink -f "$0")")"
  output_dir="$script_dir/dependencies"
fi

if [ -z "$humandb" ]; then
  # Set humandb to output_dir
  humandb="$output_dir"
fi

# Use the values of the options as needed
echo "Annovar: $annovar"
echo "Output directory: $output_dir"
echo "Annovar Human DB: $humandb"
echo "Config File: $config_file"

# Download reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P $output_dir
gunzip $output_dir/hg38.fa.gz
samtools faidx $output_dir/hg38.fa

# Download other data requirements
wget https://ftp.ncbi.nih.gov/gene/DATA/mim2gene_medgen -P $output_dir
wget https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz -P $output_dir
gunzip $output_dir/Homo_sapiens.gene_info.gz
wget http://purl.obolibrary.org/obo/hp/hp.obo -P $output_dir
wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa -P $output_dir
wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt -P $output_dir

# Download annovar databases
$annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene $humandb
$annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar gnomad30_genome $humandb
$annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar clinvar_20220320 $humandb
$annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar dbnsfp42a $humandb


# Update the YAML config file with absolute paths
if [ -f "$config_file" ]; then

  # Create a temporary file for editing
  temp_file=$(mktemp)
  
  # Use sed to make the replacements and save to the temporary file
  sed "s|annovar_scripts:.*|annovar_scripts: $annovar|" "$config_file" |
  sed "s|annovar_db:.*|annovar_db: $humandb|" |
  sed "s|reference_path:.*|reference_path: $output_dir/hg38.fa|" |
  sed "s|^mutpredindel:.*|mutpredindel: $output_dir/mutpredindel|" |
  sed "s|^muptredindel_MCR:.*|muptredindel_MCR: $output_dir/mutpred-indel-mcr|" |
  sed "s|gene_info:.*|gene_info: $output_dir/Homo_sapiens.gene_info|" |
  sed "s|hpo:.*|hpo: $output_dir/hp.obo|" |
  sed "s|hpoa:.*|hpoa: $output_dir/phenotype.hpoa|" |
  sed "s|phenotype_gene:.*|phenotype_gene: $output_dir/phenotype_to_genes.txt|" |
  sed "s|mim2gene:.*|mim2gene: $output_dir/mim2gene_medgen|" |
  sed "s|common_variants:.*|common_variants: $data_dir/common_variants.csv|" |
  sed "s|pheno_score_source:.*|pheno_score_source: $output_dir/HPO_with_gene_IC.csv|" > "$temp_file"

  # Replace the original config file with the temporary file
  mv "$temp_file" "$config_file"
  
  echo "Updated the config file: $config_file"
else
  echo "Config file not found: $config_file"
fi
