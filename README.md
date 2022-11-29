# CAVaLRi
CAVaLRi (Clinical Assessment of Variants by Likelihood Ratios) was designed to better prioritize clinically actionable variants in genetic disease. For information regarding the architecture and motivation behind CAVaLRi, please refer to our publication: {link}

## Dependencies
CAVaLRi requries a series of python packages, variant annotation, and LIRICAL dependencies in order to run properly.

### Conda
CAVaLRi utilizes a conda environment to manage Python packages. The following code will download and run the most recent installer for Miniconda. If conda is already installed on your machine, skip this step.
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
Miniconda3-latest-Linux-x86_64.sh
```

Be sure that the conda binary has been added to your path. To create the CAVaLRi conda environment, run the following command from the repository root:
```
conda env create -f dependencies/conda_environment.yml
```

Once this is created, you can activate the conda environment with the following command:

NOTE: CAVaLRi will not run to completion of the cavalri conda environment is not active

### Variant Annotation Databases
CAVaLRi utilizes the ANNOVAR variant annotation library to functionally annotate inputted genetic variants. The RefSeq, gnomAD, and clinvar databases are required. These databases can be downloaded by running the following command:
```
bash get_dependencies.sh
```

This script will download not only the ANNOVAR required databases, but all LIRICAL dependencies as well.

NOTE: These databases require ~80GB of storage

### Reference Genome
The hg38 reference genome is required to normalize inputted VCF files. This reference genome is downloaded from the UCSC servers, and is also included as part of the get_dependencies.sh bash script.

NOTE: The unzipped reference genome requries ~4GB of storage

### LIRICAL
CAVaLRi extends the LIRICAL likelihood ratio framework, and thus requires LIRICAL to run for each provided case. As such, the dependencies of LIRICAL are inherited by CAVaLRi. These depenencies are downloaded along with the variant annotation databases by the get_dependencies.sh bash script.

NOTE: These LIRICAL dependencies require ~50GB of storage

## Running CAVaLRi
Once the necessary dependencies are installed, CAVaLRi can be run via two wrapper scripts that create the necessary CAVaLRi objects to process the inputted cases. Inputted cases are represented via a JSON file with the following structure (see example/example.json).
```
{
    "phenotype":"example.pheno.tsv",
    "vcf":"example.vcf.gz",
    "biological_sex":"M",
    "proband":"example_proband",
    "mother":"example_mother",
    "father":"example_father",
}
```

The "phenotype" file is a one-column csv file with the header "HPO ID" that contains a list of all phenotypes associated with the patient (or proband).

The "vcf" file is a multi-sample vcf file containing genotype information from the proband, the mother, and the father. The "proband", "mother", and "father" input keys designate which sample refers to the proband, the mother, and the father, respectively.

NOTE: CAVaLRi currently only supports VCF files aligned to hg38. Be sure to liftover hg19 VCFs prior to running CAVaLRi.

The "biological_sex" field indicates the biological sex of the proband (either M or F)

Once these input files are prepared for each case, we are ready to run CAVaLRi

NOTE: Validation logic will be applied to the filed indicated in the input fields, and any validation errors will be raised prior to running CAVaLRi

### Run a single case
To run CAVaLRi on a single case, the cavalri_case.py script should be utilized. To learn more about the inputs to this file run:
```
python cavalri_case.py -h
```
Two arguments are available, --input and --output_dir. The --input argument is required, and should indicate the CAVaLRi json input file for the case. The --output_dir argument is optional. If provided, the output files from the case will be redirected to the indicated directory. Otherwise, the output will be saved in the same directory that houses the input file.

### Run multiple cases
To run CAVaLRi on multiple subjects at once, the cavalri.py script should be utilized. To learn more about the inputs to this file run:
```
python cavalri_cohort.py -h
```

Two arguments are available, --input_dir and --output_dir. The --input_dir argument is required, and should indicate a path where multiple CAVaLRi json input files are located. The --output_dir argument is optional. If provided, the output files from the case will be redirected to the indicated directory. Otherwise, the output will be saved in the same directory indicated by the --input_dir argument.

## Interpreting the output
CAVaLRi outputs two files per case,
1. A JSON file containing all data used in calculating the CAVaLRi post-test probability. This JSON file is keyed by disease, and contains the associated phenotype, genotype, and mode of inhertence likelihood ratios.
2. A csv summary file ranking the candidate disease genes by diagnostic post-test probability. This file is meant for variant scientists and clinicians to prioritize likely diagnostic variants.

When running cavalri_cohort.py, additional summary files are provided, including a cohort_summary.csv file. If any known diagnostic variants are provided, plots are generated to visualize CAVaLRi performance.
