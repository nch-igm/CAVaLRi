import os
import re
import pandas as pd
import json
import pickle


def build_case_data(case):

    config = case.cohort.config
    root_path = case.cohort.root_path

    # Get variants
    pathogenic_df = case.genotype.pathogenic_variants.copy()
    pathogenic_df['TYPE'] = 'SHORT'
    pathogenic_cnv_df = case.genotype.pathogenic_cnvs.copy()
    pathogenic_cnv_df['TYPE'] = 'CNV'
    pathogenic_cnv_df = pathogenic_cnv_df.rename(columns = {'CHROMOSOME':'CHROM'})
    pathogenic_cnv_df['REF'] = pathogenic_cnv_df['END']
    pathogenic_cnv_df['ALT'] = pathogenic_cnv_df['CHANGE']
    # var_df = case.genotype.pathogenic_variants.copy()

    # Intialize case data
    case_data = {'genes': {k:{} for k in 
            set(pathogenic_df['GENE_ID']).union(set(gene for sublist in pathogenic_cnv_df['INTERSECTING_GENES'] for gene in sublist))
        }}
    
    # Read in gene_disease dataframe
    gene_df = pd.read_csv(config['gene_info'], sep = '\t')
    gene_lookup = {int(row['GeneID']):row['Symbol'] for idx, row in gene_df.iterrows()}
    gene_df = gene_df[['GeneID','Symbol','Synonyms']].rename(columns = {'Symbol':'geneSymbol'}).astype({'GeneID':str})
    gene_df['Synonyms_'] = gene_df['Synonyms'].str.split('|')
    gene_syn_df = gene_df.explode('Synonyms_', ignore_index=True)[['GeneID','Synonyms_']].rename(columns = {'Synonyms_':'geneSymbol'})
    gene_syn_df = gene_syn_df[gene_syn_df['geneSymbol'] != '-'].reset_index(drop=True)
    gene_syn_df = gene_syn_df[~gene_syn_df['geneSymbol'].isin(gene_df['geneSymbol'])]
    gene_df = pd.concat([gene_df[['GeneID','geneSymbol']], gene_syn_df]).sort_values('GeneID').reset_index(drop=True).astype({'GeneID':int})
    gene_df = gene_df[~gene_df['geneSymbol'].str.startswith('LOC')].reset_index(drop=True)
    gene_df = gene_df.groupby('geneSymbol').min().reset_index()[['GeneID','geneSymbol']].astype({'GeneID':int})

    gene_disease_df = pd.read_csv(config['mim2gene'], sep = '\t').rename(columns = {'#MIM number':'OMIM'})
    gene_disease_df = gene_disease_df[(gene_disease_df['GeneID'] != '-') & (gene_disease_df['type'] == 'phenotype')].astype({'GeneID':int})
    gene_disease_df = gene_disease_df.merge(gene_df, on = 'GeneID')

    # Remove disease without an annotated mode of inheretence
    moi_df = case.cohort.moi.copy()
    moi_df.to_csv('/igm/projects/CAVaLRi/example/case/moi.csv', index = False)
    moi_diseases = [int(moi[moi.find(':')+1:]) for moi in list(set(moi_df['omimId']))]
    gene_disease_df = gene_disease_df[gene_disease_df['OMIM'].isin(moi_diseases)]

    for k,v in case_data['genes'].items():

        omim_ids = gene_disease_df[gene_disease_df['GeneID'] == k]['OMIM'].to_list()
        for oi in omim_ids:
            v[str(oi)] = {
                'MOI': ','.join(moi_df[moi_df['omimId'] == f'OMIM:{oi}']['moi'].to_list())
            }
    
        # Add genotype data
        g_short_df = pathogenic_df[pathogenic_df['GENE_ID'] == k]
        g_cnv_df = pathogenic_cnv_df[pathogenic_cnv_df['INTERSECTING_GENES'].apply(lambda x: k in x)]
        g_cnv_df['POS'] = ''
        g_cnv_df['score'] = ''

        # cols = ['CHROM','POS','REF','ALT','proband','mother','father','score','TYPE']
        g_df = pd.concat([g_short_df, g_cnv_df])
        v['gene_data'] = {
            'symbol': gene_lookup[k],
            'variants': [

                {
                    'hg38_position': f"{row['CHROM']}:{row['POS']}{row['REF']}>{row['ALT']}",
                    'score': row['score'],
                    'clinvar_pathogenic': row['clinvar_path_sig'],
                    'function': row['func'],
                    'exon_function': row['exon_func'],
                    'proband': json.loads(row['proband'])['GT'],
                    'mother': json.loads(row['mother'])['GT'] if case.mother != 'Unavailable' else 'Unavailable',
                    'father': json.loads(row['father'])['GT'] if case.father != 'Unavailable' else 'Unavailable'
                }

                if row['TYPE'] == 'SHORT' else 

                 {
                    'hg38_position': f"{row['CHROM']}:{int(row['START'])}-{int(row['END'])}_{row['CHANGE']}",
                    'score': row['PATHOGENIC'],
                    'intersecting_genes': ','.join([str(x) for x in row['INTERSECTING_GENES']]) if len(row['INTERSECTING_GENES']) > 0 else 'None',
                    'proband': json.loads(row['proband'])['GT'],
                    'mother': json.loads(row['mother'])['GT'] if case.mother != 'Unavailable' else 'Unavailable',
                    'father': json.loads(row['father'])['GT'] if case.father != 'Unavailable' else 'Unavailable'
                }
                
            for idx,row in g_df.iterrows()]
        }
    
    case_data = {'genes': {k:v for k,v in case_data['genes'].items() if len(v.keys()) > 1}}
    case_data = {'genes': {k:{k_:v_ for k_,v_ in v.items() if
                                k_ == 'gene_data'
                                    or
                                re.search('AR',v_['MOI']) and len(v['gene_data']['variants']) >= 2
                                    or
                                re.search('AD',v_['MOI']) and len(v['gene_data']['variants']) >= 1
                            } for k,v in case_data['genes'].items() # if 
                        }}
    case.case_data = case_data
    
    