import sys
import os
import re
import pandas as pd
import json
import obonet


def non_zero(freq, F_d):
    return freq if freq != 0 else min([x for x in F_d.values() if x != 0])


def get_hpoa_disease_frequency(query_term, F_d, _genes):
    """
    F_d: subset of annotated HPO frequencies in disease d
    """
    
    # First the case when the disease has an associated phenotype at least as 
    # specific as the query term -- term is (an ancestor of) a disease term:
    if query_term in F_d:
        return 'Ancestraly Closed', non_zero(F_d[query_term], F_d)

    # Now suppose query term is not in ancestral closure of disease term set;
    # in this case we'll take the lowest frequency found among the most 
    # specific common ancestors of the query term and the disease pheno set
    else:

        common_ancestors = {ca:F_d[ca] for ca in hpo.lower_boundary( hpo.ancestors(query_term) & set(F_d) )}
        best_common_ancestor = None
        best_ca_score = 0
        for ca_k, ca_v in common_ancestors.items():
            disease_freq = ca_v
            query_genes = 1 if len(_genes[query_term]) == 0 else len(_genes[query_term])
            ca_genes = 1 if len(_genes[ca_k]) == 0 else len(_genes[ca_k])
            penalty = query_genes / ca_genes
            score = disease_freq * penalty
            score = disease_freq
            if score > best_ca_score:
                best_common_ancestor = ca_k
                best_ca_score = score
        
        return best_common_ancestor, non_zero(best_ca_score, F_d)



def get_hpoa_background_frequency(pheno, bkgd_freq):
    
    #TODO Add variables that condition on the background frequency
    try:
        return bkgd_freq[pheno]

    except KeyError:
        print(f"""
            There is not a background frequency available for {pheno}. Taking
            the minimum frequency available in the provided background
            frequencies.
            """)
        return 0


def score_disease_phenotype(F_d, bkgd_freq, _genes, case):

    # Initialize result
    phenos = {p:{} for p in case.phenotype.phenotypes.keys()}
    
    # Add each term to the result
    for pheno in case.phenotype.phenotypes:
        phenos[pheno]['common_ancestor'], phenos[pheno]['disease_freq'] = get_hpoa_disease_frequency(pheno, F_d, _genes)
        if phenos[pheno]['common_ancestor'] == 'Ancestraly Closed':
            background_freq = get_hpoa_background_frequency(pheno, bkgd_freq)
        elif phenos[pheno]['common_ancestor'] == None:
            background_freq = 1
        else:
            background_freq = get_hpoa_background_frequency(phenos[pheno]['common_ancestor'], bkgd_freq)
        phenos[pheno]['background_freq'] = min([bg for bg in bkgd_freq.values() if bg > 0]) if background_freq == 0 else background_freq
        phenos[pheno]['LR'] = phenos[pheno]['disease_freq'] / phenos[pheno]['background_freq']

        # Impose a penalty if the most recent common ancestor is the root node
        if phenos[pheno]['common_ancestor'] == 'HP:0000118':
            phenos[pheno]['LR'] = case.cohort.config['pheno_root_penalty']

    return phenos


def score_phenotypes(case):

    config = case.cohort.config
    sys.path.append(config['hpo_walk_dir'])
    from hpo_walk.dag import ontology
    from hpo_walk.annotations import get_phenotype_disease_gene_df, build_propagated_frequency_map

    # Read in HPO
    global hpo
    hpo = ontology(os.path.join(case.cohort.root_path, config['hpo']))

    # Intialize HPO cohort frequencies
    # with open(os.path.join(case.cohort.root_path, config['hpo_background']), 'r') as json_file:
    with open(os.path.join(case.cohort.root_path, config['hpo_bkgd_frequencies']), 'r') as json_file:
        bkgd_freq = json.load(json_file)
    
    # Read in the HPOA
    hpoa_path = os.path.join(case.cohort.root_path, config['hpoa'])
    gene_disease_path = os.path.join(case.cohort.root_path, config['phenotype_gene'])

    # Intialize frequency map
    F = build_propagated_frequency_map(hpoa_path)

    # Read in information content scores for all phenotypes
    hpo_ic_path = os.path.join(case.cohort.root_path, config['pheno_score_source'])
    hpo_ic_df = pd.read_csv(hpo_ic_path)

    # Get gene counts for all HPO terms
    _genes = {x: set() for x in hpo.terms}
    pdg = get_phenotype_disease_gene_df(gene_disease_path)
    for x,gene in zip(pdg['hpo_id'],pdg['gene']):
        for u in hpo.ancestors(x):
            try:
                _genes[u].add(gene)
            except KeyError: 
                pass

    # Score the phenotypes for every disease
    for v in case.case_data['genes'].values():
        for omim in v.keys():
            
            # Get scores for all phenotype terms
            if f'OMIM:{omim}' in F.keys():
                v[omim]['phenotype_scores'] = score_disease_phenotype(F[f'OMIM:{omim}'], bkgd_freq, _genes, case)

    return case.case_data