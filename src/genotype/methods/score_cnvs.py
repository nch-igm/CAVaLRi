import sys
import os
import re
import json
import vcf
import pandas as pd
sys.path.append('/igm/projects/')
import CNVoyant as cnv

def read_cnvs(genotype):

    vcf_reader = vcf.Reader(filename = genotype.cnv_path, compressed=True, encoding='ISO-8859-1')

    # Initialize list to capture variants
    var_list = []

    # Set sample list
    samples = {
        'proband':genotype.case.proband,
        'mother': genotype.case.mother,
        'father': genotype.case.father
    }

    for var in vcf_reader:
        
        # Variant specific
        start_pos = 3 if re.search('chr', var.CHROM) else 0
        chrom = var.CHROM[start_pos:]
        start = var.POS

        stop = var.INFO['END']
        type_ = var.INFO['SVTYPE']

        var_row = [chrom, start, stop, type_]
        columns = ['CHROMOSOME','START','END','CHANGE']
        

        # Sample specific
        for k, sample in samples.items():
            for s in var.samples:
                if s.sample == sample:

                    # Try to determine AD, AF, and DP
                    try:
                        af = round(0 if s.data.AD[1] == 0 else float(s.data.AD[1]/s.data.DP), 3)
                    except:
                        af = None
                    
                    try:
                        ad = s.data.AD
                    except:
                        ad = None

                    try:
                        dp = s.data.DP
                    except:
                        dp = None

                    # Add sample data and add sample name to columns
                    var_row.append(json.dumps({
                        'GT': s.data.GT
                    }))
                    columns.append(k)

        # Add variant to list, which will be converted back into a data frame
        var_list.append(var_row)

    if len(var_list) > 0:
        return pd.DataFrame(var_list, columns = columns, dtype=str)
    
    return pd.DataFrame(columns  = ['CHROMOSOME','START','END','CHANGE'])


def score_cnvs(genotype):
    
    # Read in VCF records
    cnv_df = read_cnvs(genotype)
    
    # Get dependencies if necessary
    db = cnv.DependencyBuilder()
    db.build_all()

    # Get features
    if len(cnv_df.index) > 0:
        fb = cnv.FeatureBuilder(variant_df = cnv_df, 
                                data_dir = db.data_dir,
                                ref = genotype.case.cohort.config['reference_path'])
        fb.get_features()

        # Get predictions
        cl = cnv.Classifier(norm = True)
        genotype.cnvs = cl.predict(input = fb.feature_df)

    else:
        genotype.cnvs = pd.DataFrame(columns = list(cnv_df.columns) + ['PATHOGENIC','proband','INTERSECTING_GENES'] )
