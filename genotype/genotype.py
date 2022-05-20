from .read_variants import *
from .annotate_variants import *
from .filter_variants import *

class Genotype:
    """
    A genotype exists for every position within the genome,
    and this class is designed to capture the characteristics
    surrounding each variant with respect to the reference 
    genome.

    To control versioning of annotation databases, variants will
    be annotated in the NCH IGM application, Varhouse. Genotypes
    may be extracted from standardized input files (.vcf) and 
    consequently sent to an AWS EMR cluster for annotation. 

    Annotations will be added in aggregate with other child cases
    in a cohort. This annotated variant file is then referenced
    in the methods described.
    """
    
    # INITIALIZE
    def __init__(self, case, genotype_path, genome_build = 'hg38'):
        self.genome_build = genome_build
        self.genotype_path = genotype_path
        self.case = case

    # METHODS
    def read_variants(self):
        self.variants = read_variants(self.genotype_path)
    
    def annotate_variants(self):
        self.annotated_variants = annotate_variants(self)

    def filter_variants(self):
        self.filtered_variants = filter_variants(self)
