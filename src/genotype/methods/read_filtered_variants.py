import os
import re
import json
import shlex
import pandas as pd
import vcf
import subprocess


def worker(cmd):
    parsed_cmd = shlex.split(cmd)
    p = subprocess.Popen(parsed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    return out.decode() if out else err.decode()


def read_filtered_variants(genotype):

    config = genotype.case.cohort.config
    root_path = genotype.case.cohort.root_path

    # Read in vcf
    vcf_reader = vcf.Reader(filename = genotype.genotype_path, compressed=True, encoding='ISO-8859-1')

    # Initialize list to capture variants
    var_list = []

    # Set sample list
    samples = {
        'proband':genotype.case.proband,
        'mother': genotype.case.mother,
        'father': genotype.case.father
    }

    # Prepare gene data frame
    cols = ['GeneID','Symbol','Synonyms']
    gene_df = pd.read_csv(
        os.path.join(root_path, config['gene_info']),
        sep = '\t'
        )[cols]

    def re_match(row):
        return True if row['GeneID'] in row['Synonyms'].split('|') else False
    
    def get_gene_id(symbol, gene_df):
        gene_match_df = gene_df[gene_df['Symbol'] == symbol]
        l = len(gene_match_df.index)
        if l == 1:
            return str(int(gene_match_df.reset_index(drop=True).loc[0,'GeneID']))
        elif l == 0:
            gm = gene_df[gene_df.apply(re_match, axis = 1)]
            if len(gm.index) > 0:
                return str(int(gm.reset_index(drop=True).loc[0,'GeneID']))
            else:
                return None
        else:
            return str(int(gene_match_df.reset_index(drop=True).loc[0,'GeneID']))

    for var in vcf_reader:
        
        # Variant specific
        start_pos = 3 if re.search('chr', var.CHROM) else 0
        chrom = var.CHROM[start_pos:]
        pos = var.POS
        ref = var.REF
        alt = ','.join([str(i) for i in var.ALT])
        info = json.dumps(var.INFO)

        # Get NCBI Gene ID
        gene = var.INFO['Gene.refGene'][0]
        gene = gene if not re.search('&', gene) else gene[:gene.find('&')]
        ncbi_id = get_gene_id(gene, gene_df)
        var_row = [chrom, pos, ref, alt, gene, ncbi_id, info]
        columns = ['CHROM', 'POS', 'REF', 'ALT', 'GENE', 'GENE_ID', 'INFO']
        

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
                        'GT': s.data.GT,
                        'AD': ad,
                        'AF': af,
                        'DP': dp
                    }))
                    columns.append(k)

        # Add variant to list, which will be converted back into a data frame
        var_list.append(var_row)
    
    df = pd.DataFrame(var_list, columns = columns, dtype=str)
    df = df[df['GENE_ID'].notna()]

    return df
