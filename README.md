# CAVaLRi
CAVaLRi (Clinical Assessment of Variants by Likelihood Ratios) was designed to better prioritize potential diagnostic genes in rare genetic diseases. For information regarding the architecture and motivation behind CAVaLRi, please refer to our publication: {link}. The figure below illustrates the additions to the original likelihood ratio framework defined in LIRICAL (Robinson et al.). These extensions include (1) filtering phenotypes to only consider the most genetically informative, (2) incorporating parental genotypes with a novel segregation likelihood ratio, and (3) weighting the likelihood ratios by relative importance.
![image](https://github.com/nch-igm/CAVaLRi/assets/72405035/999b82e8-ac8a-4d96-826a-d24e3b9e6b9a)

## Installation Introduction for CAVaLRi

1. *Navigate to the Desired Directory*

First, we recommend creating and moving to an apps directory within your home folder. This will be the location where the CAVaLRi software will reside. To do this, open your terminal and run:
```
mkdir -p ~/apps
cd ~/apps
```

2. *Clone the CAVaLRi Repository*

Before you begin the installation process, ensure you have git installed on your system. With the terminal now pointed to the ~/apps directory, you can clone the CAVaLRi repository:
```
git clone https://github.com/nch-igm/CAVaLRi.git
```

3. *Move into new directory*

At this point, you’ve successfully cloned the CAVaLRi repository to your system. Navigate to the CAVaLRi directory to continue installing the software:
```
cd CAVaLRi
```

## Prerequisites

### Git
Ensure you have git installed on your system. If not, you can install it using the package manager of your operating system.

### Mamba
CAVaLRi utilizes a conda environment to manage Python packages. Conda takes quite a long time to install dependencies, so we recommend utilizing mamba instead (micromamba for minimal installation). The following code will download and run the most recent installer for micromamba [see micromamba configuration instructions](https://mamba.readthedocs.io/en/latest/micromamba-installation.html).
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Once configured, create the CAVaLRi conda environment by running the following command from the repository root:

Note: Mac users with an M1 processor should utilize the M1 specific environment file (dependencies/mamba_environment_M1.yml)
```
micromamba env create -f dependencies/mamba_environment.yml
```

Once this is created, you can activate the conda environment with the following command:
```
micromamba activate cavalri
```

NOTE: CAVaLRi will not run to completion if the cavalri conda environment is not active

### Annovar (Variant Annotation)
CAVaLRi utilizes the ANNOVAR variant annotation library to functionally annotate inputted genetic variants. ANNOVAR must be requested and downloaded manually, see instructions [here](https://annovar.openbioinformatics.org/en/latest/user-guide/download/).

## Dependencies
CAVaLRi requires a reference genome, variant annotation databases, and phenotype-disease-gene associations in order to run properly.

### Reference Genome
The hg38 reference genome is required to normalize inputted VCF files. This reference genome is downloaded from the UCSC servers and is also included as part of the get_dependencies.sh bash script.

### Variant Annotation Databases
The RefSeq, gnomAD, ClinVar, and dbNSFP databases are required. These databases can be downloaded by running the following command:
```
bash get_dependencies.sh -a /path/to/annovar
```

This script will download ANNOVAR required databases, disease-gene-phenotype annotations, and ontology files.

NOTE: These dependencies together require ~120GB of storage. Also, the database versions in get_dependencies.sh can be updated manually to specific versions if needed.

## Running CAVaLRi
Once the necessary dependencies are installed, CAVaLRi can be run via two wrapper scripts that create the necessary CAVaLRi objects to process the inputted cases. Inputted cases are represented via a JSON file with the following structure (see example/example.json).
```
{
    "phenotype": "example/case/example.pheno.csv",
    "vcf": "example/case/example.vcf.gz",
    "pedigree": "example/case/example.ped",
    "proband": "PROBAND"
}
```

The “phenotype” argument indicates a one-column csv file with the header “HPO ID” that contains a list of all phenotypes associated with the patient (or proband). Phenotypes must be represented using HPO IDs, e.g., HP:0001252, with each new phenotype on a new line in the file. 

The "pedigree" argument indicates a .ped file that designates the mother and father samples as well as the affected statuses of each individual in the pedigree. For more details, see [.PED formatting](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).

The "vcf" argument indicates a multi-sample vcf file containing genotype information from the proband, the mother, and the father. The individual identifiers designating which sample refers to the proband, the mother, and the father, should be identified in the provided .ped file.

NOTE: CAVaLRi currently only supports VCF files aligned to hg38. Be sure to liftover hg19 VCFs prior to running CAVaLRi.

The "proband" argument indicates the name of the proband sample in the vcf


Once these input files are prepared for each case, we are ready to run CAVaLRi

NOTE: Validation logic will be applied to the data indicated in the input fields, and any validation errors will be raised prior to running CAVaLRi

### Run a single case
To run CAVaLRi on a single case, the cavalri_case.py script should be utilized.
```
python cavalri_case.py -i example/case/example.json -o testing
```
Two arguments are available, --input and --output_dir. The --input argument is required and should indicate the CAVaLRi json input file for the case. The --output_dir argument is optional. If provided, the output files from the case will be redirected to the indicated directory. Otherwise, the output will be saved in the same directory that houses the input file.

### Run multiple cases
To run CAVaLRi on multiple subjects at once, the cavalri.py script should be utilized.
```
python cavalri_cohort.py -i example/cohort -o testing
```

Two arguments are available, --input_dir and --output_dir. The --input_dir argument is required and should indicate a path where multiple CAVaLRi json input files are located. The --output_dir argument is optional. If provided, the output files from the case will be redirected to the indicated directory. Otherwise, the output will be saved in the same directory indicated by the --input_dir argument.

## Interpreting the output
CAVaLRi outputs two files per case,
1. A JSON file containing all data used in calculating the CAVaLRi posterior probability. This JSON file contains two primary keys, "phenotypes" and "genes". A list of the patient's scored phenotypes is stored in the "phenotype" values. The "genes" value contains all genotype, phenotype, and segregation data needed to calculate the CAVaLRi score. This value is itself a dictionary keyed by gene (NCBI ID) whose values are representations of OMIM diseases associated with the keyed gene. This file is provided if more granular analysis is required.
2. A csv summary file ranking the candidate disease genes by diagnostic posterior probability. This file is meant for variant scientists and clinicians to prioritize likely diagnostic variants. Of note, likelihood ratio values are log-transformed for interpretability. The values displayed in genoLR and segLR are scaled according to the optimal exponents recorded in the configuration file.

When running cavalri_cohort.py, an additional summary file is generated that sorts all candidate genes in the cohort by diagnostic posterior probability.

## References
1. Robinson PN, Haendel MA. Ontologies, Knowledge Representation, and Machine Learning for Translational Research: Recent Contributions. Yearb Med Inform. 2020;29:159–62.
