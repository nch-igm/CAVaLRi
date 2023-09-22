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


def get_gene_id(symbol, gene_df): 
    try:
        return gene_df[gene_df['geneSymbol'] == symbol].reset_index().loc[0,'GeneID']
    except:
        return False 

def get_gene_symbol(gene_id, gene_df_):
    try:
        return gene_df_[gene_df_['GeneID'] == gene_id].reset_index().loc[0,'geneSymbol']
    except:
        return False   

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

    # Prepare gene data frames
    gene_df = pd.read_csv(os.path.join(root_path, config['gene_info']), sep = '\t')
    gene_df = gene_df[['GeneID','Symbol','Synonyms']].rename(columns = {'Symbol':'geneSymbol'}).astype({'GeneID':str})
    gene_df_ = gene_df.copy()[['GeneID','geneSymbol']].astype({'GeneID':int})
    gene_df['Synonyms_'] = gene_df['Synonyms'].str.split('|')
    gene_syn_df = gene_df.explode('Synonyms_', ignore_index=True)[['GeneID','Synonyms_']].rename(columns = {'Synonyms_':'geneSymbol'})
    gene_syn_df = gene_syn_df[gene_syn_df['geneSymbol'] != '-'].reset_index(drop=True)
    gene_syn_df = gene_syn_df[~gene_syn_df['geneSymbol'].isin(gene_df['geneSymbol'])]
    gene_df = pd.concat([gene_df[['GeneID','geneSymbol']], gene_syn_df]).sort_values('GeneID').reset_index(drop=True).astype({'GeneID':int})
    gene_df = gene_df[~gene_df['geneSymbol'].str.startswith('LOC')].reset_index(drop=True)
    gene_df = gene_df.groupby('geneSymbol').min().reset_index()


    # def re_match(row, symbol):
    #     return True if symbol in row['Synonyms'].split('|') else False
    
    # def get_gene_id(symbol, gene_df):
    #     gene_match_df = gene_df[gene_df['Symbol'] == symbol]
    #     l = len(gene_match_df.index)
    #     if l == 1:
    #         return str(int(gene_match_df.reset_index(drop=True).loc[0,'GeneID']))
    #     elif l == 0:
    #         gm = gene_df[gene_df.apply(re_match, symbol = symbol, axis = 1)]
    #         if len(gm.index) > 0:
    #             return str(int(gm.reset_index(drop=True).loc[0,'GeneID']))
    #         else:
    #             return None
    #     else:
    #         return str(int(gene_match_df.reset_index(drop=True).loc[0,'GeneID']))

    # Get OMIM disease gene IDs
    omim_gene_ids = pd.read_csv(os.path.join(root_path, config['mim2gene']), sep = '\t')
    omim_gene_ids = omim_gene_ids[omim_gene_ids['GeneID'] != '-'].astype({'GeneID':int})
    omim_gene_ids = omim_gene_ids.merge(gene_df, on = 'GeneID')
    omim_gene_ids = list(set(omim_gene_ids[omim_gene_ids['type'] == 'phenotype']['GeneID']))

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

        # One annotated gene
        if not re.search('x3b', gene):
            gene_id = get_gene_id(gene, gene_df)

        # More than one annotated gene
        else:
            gl_ = []
            gl = gene.split('\\x3b') # Gene list
            for gs in gl:
                gl_id = get_gene_id(gs,gene_df)
                if gl_id:
                    gl_.append(gl_id)
            
            # Break if there is no match
            if len(gl_) == 0:
                continue

            gene_id = gl_[0]

            aa = var.INFO['AAChange.refGene'][0] # First amino acid change listed
            if aa:
                aa_pos = aa.find(':')
                if aa_pos != -1:
                    aa = aa[:aa_pos]
                    aa_id = get_gene_id(aa, gene_df)
                    if aa_id:
                        print(aa_id, gl_, True if aa_id in omim_gene_ids else False)
                        if aa_id in gl_ and aa_id in omim_gene_ids:
                            gene_id = aa_id

            else:
                omim_ = [g for g in gl_ if g in omim_gene_ids] # OMIM genes
                gene_id = omim_[0] if len(omim_) > 0 else gl_[0]


        if gene_id not in gene_df['GeneID'].to_list():
            continue

        # ncbi_id = gene_df[gene_df['GeneID'] == gene_id].reset_index().loc[0,'GeneID']
        gene_symbol = get_gene_symbol(gene_id, gene_df_)
        var_row = [chrom, pos, ref, alt, gene_symbol, gene_id, info]
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
    df = df[df['GENE_ID'].notna()].astype({'GENE_ID': int})

    return df
