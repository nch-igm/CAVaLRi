import sys
import os
import re
import pandas as pd
import json

sys.path.append('..')
# from hpo_walk import ontology
# from hpo_walk.annotations import get_phenotype_disease_gene_df, build_propagated_frequency_map
from config import *

# Intialize HPO frequency map
global hpo_map
hpo_map = pd.read_csv(os.path.join(config['project_root'], config['hpo_frequency_map']))

# Intialize HPO cohort frequencies
global bkgd_freq
with open(os.path.join(config['project_root'], config['hpo_background']), 'r') as json_file:
    bkgd_freq = json.load(json_file)

# Intialize frequency map
global F
F = build_propagated_frequency_map()


def hpo_distance(hpo, term1, term2):
    max_hpo_length = 20
    counter, frontier = 0, {term1}
    while term2 not in frontier and counter < max_hpo_length:
        frontier = hpo.neighbors(frontier)
        counter += 1
    return counter if counter != max_hpo_length else 999


def parse_frequency(freq):

    # If the frequency is null
    if pd.isna(freq):
        return 0.0999
    
    # If the frequency is an HPO term
    if re.search('HP:', freq):
        if config['hpo_frequency_bound'] == 'upper':
            return hpo_map.loc[hpo_map['HPO_ID'] == freq].reset_index(drop=True)['upper_bound'][0]
        else:
            return hpo_map.loc[hpo_map['HPO_ID'] == freq].reset_index(drop=True)['lower_bound'][0]
    
    # If the frequency is a fraction
    if re.search('/', freq):
        frac = freq.split('/')
        return int(frac[0])/int(frac[1])
    
    # If the frequency is a percentage
    if re.search('%', freq):
        return int(freq[:-1])/100


def get_hpoa_disease_frequency(query_term, F_d, hpo):


# def get_phenotype_disease_frequency(term, disease, hpo = HPO, freq_map = F):
    
    ### first the case when the disease has an associated phenotype at least as specific as the query term -- term is (an ancestor of) a disease term:
    if query_term in F_d:
        return F_d[query_term]

    ### now suppose query term is not in ancestral closure of disease term set;
    ### in this case we'll take the lowest frequency found among the most specific common ancestors of the query term and the disease pheno set
    else:
        return min(
            [F_d[u] for u in hpo.lower_boundary( hpo.ancestors(query_term) & set(F_d) )]
            )

    # Define HPO disease terms
    disease_terms = disease_term_df['HPO_ID'].unique()
    
    # Check to see if the term is directly annotated in 
    direct_match_df = disease_term_df.loc[disease_term_df['HPO_ID'] == query_term]
    if len(direct_match_df.index) == 1:
        frequency = direct_match_df.reset_index(drop=True)['Frequency'][0]
        return parse_frequency(frequency)
    # else:
    #     # Non-match scenario
    #     some_value = 1/100
    #     return some_value

    # Check for ascendents
    ascendents = {}
    for disease_term in disease_terms:
        if query_term in hpo.ancestors(disease_term):
            freq = disease_term_df.loc[disease_term_df['HPO_ID'] == disease_term].reset_index(drop=True)['Frequency'][0]
            ascendents[disease_term] = parse_frequency(freq)
    if ascendents != {}:
        return max(ascendents.values())

    # Check for descendents
    descendents = {}
    for disease_term in disease_terms:
        if query_term in hpo.descendants(disease_term):
            freq = disease_term_df.loc[disease_term_df['HPO_ID'] == disease_term].reset_index(drop=True)['Frequency'][0]
            descendents[disease_term] = parse_frequency(freq)
    if descendents != {}:
        return max(descendents.values())/2

    # Check for common ancestors
    best_common_ancestor = {
        'query_term': term,
        'disease_term': None,
        'ca_term': None,
        'frequency': 0
    }
    
    # Intialize length counter (arbitrarily large value)
    l = 1000

    # Determine common ancestors, iterate through to find the closest
    for disease_term in disease_terms:

        # Get common ancestors, remove root nodes
        common_ancestors = hpo.ancestors(term).intersection(hpo.ancestors(disease_term))
        common_ancestors.remove('HP:0000001')
        if 'HP:0000118' in common_ancestors:
            common_ancestors.remove('HP:0000118')

        for ca in common_ancestors:

            # Get frequency of disease term
            disease_freq = disease_term_df.loc[disease_term_df['HPO_ID'] == disease_term].reset_index(drop=True)['Frequency'][0]
            
            # If common ancestor is closer than currently closest ancestor
            if hpo_distance(hpo, term, ca) < l:
                best_common_ancestor = {
                    'query_term': term,
                    'disease_term': disease_term,
                    'ca_term': ca,
                    'frequency': parse_frequency(disease_freq)
                }
                l = hpo_distance(hpo, term, ca)
                
            # If common ancestor is as close as currently closest ancestor, but frequency is higher
            elif hpo_distance(hpo, term, ca) == l and parse_frequency(disease_freq) > best_common_ancestor['frequency']:
                best_common_ancestor = {
                    'query_term': term,
                    'disease_term': disease_term,
                    'ca_term': ca,
                    'frequency': parse_frequency(disease_freq)
                }
    
    if best_common_ancestor['ca_term'] is not None:
        return best_common_ancestor['frequency']


    return 0.01
                

def get_hpoa_background_frequency(term):
    
    #TODO Add variables that condition on the background frequency
    return bkgd_freq[term]


def score_phenotype(omim_id, hpo, hpoa, case):

    # Initialize result
    phenos = {}
    for p in case.phenotype.phenotypes.keys():
        phenos[p] = {}

    # Get the annotations for this disease
    disease_term_df = hpoa.loc[hpoa['DatabaseID'] == ('OMIM:' + str(omim_id))][['HPO_ID', 'Frequency']]

    # Add each term to the result
    for proband_term in case.phenotype.phenotypes:
        # phenos[proband_term]['disease_freq'] = get_hpoa_disease_frequency(proband_term, disease_term_df, hpo)
        phenos[proband_term]['disease_freq'] = get_hpoa_disease_frequency(proband_term, F[omim_id], hpo)
        phenos[proband_term]['background_freq'] = get_hpoa_background_frequency(proband_term)
        phenos[proband_term]['LR'] = max(0.01,
                phenos[proband_term]['disease_freq'] / phenos[proband_term]['background_freq'])

    return phenos


def score_phenotypes(case):

    # Read in HPO
    hpo = ontology('HPO')

    # Read in the HPOA columns
    hpoa_columns = list(pd.read_csv(
        os.path.join(config['project_root'], config['hpoa']),
        skiprows = 4,
        sep = '\t',
        nrows = 1
    ))

    # Read in the HPOA
    hpoa = pd.read_csv(
        os.path.join(config['project_root'], config['hpoa']),
        skiprows = 4,
        sep = '\t',
        usecols = [c for c in hpoa_columns if c not in ['Onset', 'Sex', 'Modifier']]
    ).rename(columns = {'#DatabaseID': 'DatabaseID'})

    # Score the phenotypes for every disease
    for gene in case.case_data['genes']:
        for d in case.case_data['genes'][gene]['diseases']:
            d['phenotype_scores'] = score_phenotype(d['omim'], hpo, hpoa, case)

    return case.case_data
