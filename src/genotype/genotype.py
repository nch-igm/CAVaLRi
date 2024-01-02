from .methods import *
import os
import time


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
    def __init__(self, case, vcf_path):
        # self.genome_build = genome_build
        # self.vcf_path = os.path.join(case.cohort.root_path, vcf_path)
        self.vcf_path = vcf_path
        self.case = case

    # METHODS
    def split_variants(self):
        split_variants(self)

    def normalize_variants(self):
        normalize_variants(self)

    def annotate_variants(self):
        annotate_variants(self)

    def filter_variants(self):
        filter_variants(self)
    
    def read_filtered_variants(self):
        read_filtered_variants(self)
    
    def score_pathogenicity(self):
        run_spliceai(self)
        run_mutpredindel(self)
        get_snv_scores(self)

    def score_cnvs(self):
        score_cnvs(self)

    def filter_pathogenic_variants(self, is_empty):
        filter_pathogenic_variants(self, is_empty = is_empty)
    
    def calculate_genoLR(self):
        calculate_genoLR(self)
