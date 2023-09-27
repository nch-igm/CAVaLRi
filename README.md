# CAVaLRi
CAVaLRi (Clinical Assessment of Variants by Likelihood Ratios) was designed to better prioritize potential diagnostic genes in rare genetic diseases. For information regarding the architecture and motivation behind CAVaLRi, please refer to our publication: {link}. The figure below illustrates the additions made to the original likelihood ratio framework defined in LIRICAL (Robinson et al.). These extensions include (1) filtering phenotypes to only consider the most genetically informative, (2) incorporating parental genotypes with a novel mode of inheritence likelihood ratio, and (3) weighting the likelihood ratios by relative importance.
![image](https://github.com/nch-igm/CAVaLRi/assets/72405035/999b82e8-ac8a-4d96-826a-d24e3b9e6b9a)


## Dependencies
CAVaLRi requries a series of python packages, variant annotation, and LIRICAL dependencies in order to run properly.

### Conda (or Mamba)
CAVaLRi utilizes a conda environment to manage Python packages. Conda takes quite a long time to install dependencies, so we recommend utilizng mamba instead (micromamba for minimal installation). The following code will download and run the most recent installer for micromamba [see micromamba configuration instructions](https://mamba.readthedocs.io/en/latest/micromamba-installation.html).
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Once configured, create the CAVaLRi conda environment by running the following command from the repository root:
```
micromamba env create -f dependencies/conda_environment.yml
```

Once this is created, you can activate the conda environment with the following command:
```
micromamba activate cavalri
```

NOTE: CAVaLRi will not run to completion if the cavalri conda environment is not active

### Variant Annotation Databases
CAVaLRi utilizes the ANNOVAR variant annotation library to functionally annotate inputted genetic variants. ANNOVAR must be requested and downloaded manually, see instructions [here](https://annovar.openbioinformatics.org/en/latest/user-guide/download/).

The RefSeq, gnomAD, ClinVar, and dbNSFP databases are required. These databases can be downloaded by running the following command:
```
bash get_dependencies.sh -a /path/to/annovar
```

This script will download not only the ANNOVAR required databases, but all other dependencies as well.

NOTE: These databases require ~190GB of storage

### Reference Genome
The hg38 reference genome is required to normalize inputted VCF files. This reference genome is downloaded from the UCSC servers, and is also included as part of the get_dependencies.sh bash script.

NOTE: The unzipped reference genome requries ~4GB of storage

## Running CAVaLRi
Once the necessary dependencies are installed, CAVaLRi can be run via two wrapper scripts that create the necessary CAVaLRi objects to process the inputted cases. Inputted cases are represented via a JSON file with the following structure (see example/example.json).
```
{
    "phenotype": "example/case/example.pheno.csv",
    "pedigree": "example/case/example.ped",
    "vcf": "example/case/example.vcf.gz",
    "proband": "PROBAND"
}
```

The "phenotype" argument indicates a one-column csv file with the header "HPO ID" that contains a list of all phenotypes associated with the patient (or proband).

The "pedigree" argument indicates a .ped file that designates the mother and father samples as well as the affected statuses of each individual in the pedgree. For more details, see [.PED formatting](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).

The "vcf" argument indicates a multi-sample vcf file containing genotype information from the proband, the mother, and the father. The individual identifiers designating which sample refers to the proband, the mother, and the father, should be identified in the provided .ped file.

NOTE: CAVaLRi currently only supports VCF files aligned to hg38. Be sure to liftover hg19 VCFs prior to running CAVaLRi.

The "proband" argument indicates the name of the proband sample in the vcf


Once these input files are prepared for each case, we are ready to run CAVaLRi

NOTE: Validation logic will be applied to the data indicated in the input fields, and any validation errors will be raised prior to running CAVaLRi

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

Three arguments are available, --input_dir, --output_dir, and --diagnostic_data. The --input_dir argument is required, and should indicate a path where multiple CAVaLRi json input files are located. The --output_dir argument is optional. If provided, the output files from the case will be redirected to the indicated directory. Otherwise, the output will be saved in the same directory indicated by the --input_dir argument. The --diagnostic_data is optional, and is meant to provide diagnostic data in cases where the diagnostic gene is known. This file is to be formatted as a csv file with two columns with headers 'CASE' and 'DIAGNOSTIC_GENE'. This diagnostic data will be used to build accuracy metrics and plots at the end of the cohort pipeline.

## Interpreting the output
CAVaLRi outputs two files per case,
1. A JSON file containing all data used in calculating the CAVaLRi post-test probability. This JSON file is keyed by disease, and contains the associated phenotype, genotype, and mode of inhertence likelihood ratios.
2. A csv summary file ranking the candidate disease genes by diagnostic post-test probability. This file is meant for variant scientists and clinicians to prioritize likely diagnostic variants.

When running cavalri_cohort.py, additional summary files are provided, including a cohort_summary.csv file. If any known diagnostic variants are provided, plots are generated to visualize CAVaLRi performance.

## Examples
To get familiarized with CAVaLRi case output
```
python cavalri_case.py --input example/case/example.json
```

To get familiarized with CAVaLRi cohort output
```
python cavalri_cohort.py --input_dir example/cohort/ --diagnostic_data example/cohort/example_diagnostic.csv
```
