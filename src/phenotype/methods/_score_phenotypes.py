import sys
import os
import re
import pandas as pd
import json
import obonet


sys.path.append('/Users/rsrxs003/projects/LIRICAL')
# from hpo_walk import ontology
from hpo_walk.dag import ontology
from hpo_walk.annotations import get_phenotype_disease_gene_df, build_propagated_frequency_map


def hpo_distance(hpo, term1, term2):
    max_hpo_length = 20
    counter, frontier = 0, {term1}
    while term2 not in frontier and counter < max_hpo_length:
        frontier = hpo.neighbors(frontier)
        counter += 1
    return counter if counter != max_hpo_length else 999


def parse_frequency(freq, config):

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


def get_hpoa_disease_frequency(query_term, F_d):
    """
    F_d: subset of annotated HPO frequencies in disease d
    """

    # def get_phenotype_disease_frequency(term, disease, hpo = HPO, freq_map = F):
    
    # First the case when the disease has an associated phenotype at least as 
    # specific as the query term -- term is (an ancestor of) a disease term:
    if query_term in F_d:
        return F_d[query_term]

    # Now suppose query term is not in ancestral closure of disease term set;
    # in this case we'll take the lowest frequency found among the most 
    # specific common ancestors of the query term and the disease pheno set
    else:
        return min(
            [F_d[ca] for ca in hpo.lower_boundary( hpo.ancestors(query_term) & set(F_d) )]
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
                

def get_hpoa_background_frequency(pheno, bkgd_freq):
    
    #TODO Add variables that condition on the background frequency
    return bkgd_freq[pheno]


def score_disease_phenotype(F_d, bkgd_freq, case):

    # Initialize result
    default_value = 1 #?
    phenos = {p:{} for p in case.phenotype.phenotypes.keys()}

    # Add each term to the result
    for pheno in case.phenotype.phenotypes:
        # phenos[pheno]['disease_freq'] = get_hpoa_disease_frequency(pheno, disease_hpoa_df, hpo)
        phenos[pheno]['disease_freq'] = get_hpoa_disease_frequency(pheno, F_d)
        phenos[pheno]['background_freq'] = get_hpoa_background_frequency(pheno, bkgd_freq)

        pLR = 1 if phenos[pheno]['background_freq'] == 0 else phenos[pheno]['disease_freq'] / phenos[pheno]['background_freq']

        phenos[pheno]['LR'] = max(0.01, pLR)

    return phenos


def score_phenotypes(case):

    config = case.cohort.config

    # Read in HPO
    global hpo
    hpo = ontology(config['hpo'])

    # Intialize HPO frequency map
    # hpo_map = pd.read_csv(os.path.join(case.cohort.root_path, config['hpo_frequency_map']))
    # hpo_map = pd.read_csv('/Users/rsrxs003/projects/LIRICAL/data/hpo_frequency_map.csv')

    # Intialize HPO cohort frequencies
    # with open(os.path.join(case.cohort.root_path, config['hpo_background']), 'r') as json_file:
    with open(config['hpo_bkgd_frequencies'], 'r') as json_file:
        bkgd_freq = json.load(json_file)

    # Intialize frequency map
    F = build_propagated_frequency_map()  #TODO Find a way to provide HPOA

    # Read in the HPOA
    # hpoa_path = os.path.join(case.cohort.root_path, config['hpoa'])
    # hpoa_path = '/Users/rsrxs003/projects/LIRICAL/data/phenotype.hpoa'
    # hpoa = pd.read_csv(hpoa_path, skiprows = 4, sep = '\t')
    # hpoa_cols = hpoa.columns
    # hpoa = hpoa[hpoa_cols].rename(columns = {'#DatabaseID': 'DatabaseID'})

    # Score the phenotypes for every disease
    for v in case.case_data['genes'].values():
        for omim in v.keys():

            # Limit HPOA to omim related terms
            # disease_hpoa_df = hpoa.loc[hpoa['DatabaseID'] == ('OMIM:' + str(omim))][['HPO_ID', 'Frequency']]
            
            # Get scores for all phenotype terms
            # F[D] needs to go into score_phenotypes
            try:
                v[omim]['phenotype_scores'] = score_disease_phenotype(F[f'OMIM:{omim}'], bkgd_freq, case)
            except Exception as e:
                # print(f'Checkpoint: 9 ERROR:{e} {time.strftime("%H:%M:%S", time.localtime())}', file = open(fo, 'a'))
                pass

    return case.case_data
