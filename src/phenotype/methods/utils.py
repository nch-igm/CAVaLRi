import sys
import os
import re
import pandas as pd
import obonet as ob
import datetime
from tqdm import tqdm
from collections import defaultdict
import networkx as nx
sys.path.append(os.path.join(__file__,'../../..'))
from config import *


# DIR = os.path.dirname(os.path.abspath(__file__))

ONT_ROOT = {
    'HPO':'HP:0000118',
    'MONDO':'MONDO:0000001',
    'ORDO':'ORPHA:C001'
    }


def get_phenotype_disease_gene_df(
    resource_dir
    # resource_dir = RESOURCE_DIR # where to load text files containing annotations
    ):
    pdg_df = pd.read_csv(
        # os.path.join(resource_dir,'phenotype_disease_gene.txt'),
        resource_dir,
        delimiter='\t',
        skiprows = 1,
        header = None)
    pdg_df.columns = [
        'hpo_id',
        'hpo_name',
        'ncbi_gene_id',
        'gene']
    return pdg_df


def normalize_frequency(freq, hpo_):
    if type(freq) == float:
        return freq
    if 'HP:' in freq: # When frequency is defined by an HPO term
        # Find the upper end of the range defined by the HPO frequency
        return max(map(int,re.findall('[\d]+',list(hpo_.synonyms(freq))[0])))/100
    if '%' in freq:
        return float(freq.strip('%'))/100
    if '/' in freq:
        x,y = freq.split('/')
        y = max(1,int(y))
        return int(x)/int(y)

def get_closure(X,directional_relative_dict):
    closure = set()
    if type(X) == str:
        closure.add(X)
        try: rels = set(directional_relative_dict[X])
        except KeyError: rels = set()
        while rels:
            r = rels.pop()
            closure.add(r)
            try: rels |= directional_relative_dict[r]
            except KeyError: pass
        return closure
    else:
        for x in X:
            closure = closure | get_closure(x, directional_relative_dict)
        return closure

def get_ontology(hpo_path, root = None):
    G =  ob.read_obo(hpo_path)
    if root:
        G = G.subgraph({root} | nx.descendants(G,root))
        return G.nodes(data = True)
    else:
        return G.nodes(data = True)

def get_dag_info(graph, root = None):
    #if root:
    #    graph = graph.subgraph({root} | nx.descendants(graph,root))
    #graph_data = graph.nodes(data=True)

    hd = dict()
    for x,x_data in graph:
        hd[x] = dict()
        hd[x]['parents'] = set()
        hd[x]['children'] = set()
        hd[x]['synonyms'] = set()
        hd[x]['xref'] = set()
        hd[x]['creation_date'] = None

    for x,x_data in graph:#_data:
        hd[x]['id'] = x
        hd[x]['name'] = x_data.get('name','')
        for y in x_data.get('synonym',set()):
            try: hd[x]['synonyms'].add(y.split('"')[1])
            except IndexError: pass
        
        hd[x]['xref'] = set(x_data.get('xref',[]))
        
        for y in x_data.get('is_a',set()):
            if y in hd:
                hd[x]['parents'].add(y)
                hd[y]['children'].add(x)
        
        hd[x]['creation_date'] = x_data.get('creation_date',None)
        for prop in x_data.get('property_value',[]):            
            if 'http://purl.org.dc.elements/1.1/date' in prop:
                hd[x]['creation_date'] = prop.split(' ')[1]
        if hd[x]['creation_date']:
            hd[x]['creation_date'] = datetime.datetime.strptime(hd[x]['creation_date'],'%Y-%m-%dT%H:%M:%SZ').date()
            


    parent_dict = {x:hd[x]['parents'] for x in hd}
    child_dict = {x:hd[x]['children'] for x in hd}
    for x in hd:
        hd[x]['ancestors'] = get_closure(x,parent_dict)
        hd[x]['descendants'] = get_closure(x,child_dict)
    
    return hd


class ontology():
    def __init__(
        self,
        hpo_path,
        source = 'HPO',
        root = None
        ):
    
        
        if (source in ONT_ROOT) and (root==None): root = ONT_ROOT[source]
        


        self.root = root
        self.data = get_dag_info(graph = get_ontology(hpo_path), root = root)
        self.terms = list(self.data)
        self._name_dict = {x: self.data[x]['name'] for x in self.terms}
        self._synonyms_dict = {x: self.data[x]['synonyms'] for x in self.terms}
        self._parents_dict = {x: self.data[x]['parents'] for x in self.terms}
        self._children_dict = {x: self.data[x]['children'] for x in self.terms}
        self._ancestors_dict = {x: self.data[x]['ancestors'] for x in self.terms}
        self._descendants_dict = {x: self.data[x]['descendants'] for x in self.terms}
        self._creation_date_dict = {x: self.data[x]['creation_date'] for x in self.terms}

    # Information retrieval methods

    def name(self,X):
        if type(X) == str:
            try:
                return self.data[X]['name']
            except KeyError:
                return X
        else:
            return {x: self.name(x) for x in X}
    
    def synonyms(self,x):
        try:
            return self.data[x]['synonyms']
        except KeyError:
            return {self.name(x)}
    
    def all_names(self,x):
        return {self.name(x)} | self.synonyms(x)
    
    def xref(self,x):
        try: return self.data[x]['xref']
        except KeyError: return set()
        
    def creation_date(self,x):
        try: return self.data[x]['creation_date']
        except KeyError: return None
    
    # Relative retrieval methods

    def parents(self,X):
        if type(X) == str:
            try: return self.data[X]['parents']
            except KeyError: return set()
        else:
            res = set()
            for x in X:
                res = res | self.parents(x)
            return res

    
    def children(self,X):
        if type(X) == str:
            try: return self.data[X]['children']
            except KeyError: return set()
        else:
            res = set()
            for x in X:
                res = res | self.children(x)
            return res
    
    def ancestors(self,X):
        return get_closure(X,self._parents_dict)
    
    def descendants(self,X):
        return get_closure(X,self._children_dict)
    
    def neighbors(self,X,n=1):
        if type(X) == str:
            return self.neighbors({X},n=n)
        if n <= 0:
            return set(X)
        if n > 0:
            start = self.neighbors(X,n-1)
            res = set(start)
            for x in start:
                for u in self.parents(x) | self.children(x):
                    res.add(u)
            return res
    
    def graph_distance(self,x,y):
        max_radius = 25  #selected to be safely larger than diameter of HPO
        dist = 0
        ball = {x}
        while y not in ball and dist <= max_radius:
            ball = self.neighbors(ball)
            dist += 1
        return dist
    
    def jaccard_desc_score(self,X,Y):
        A = self.descendants(X)
        B = self.descendants(Y)
        return len(A&B)/len(A|B)
            
    
    # Boundary/associated operations

    def leaves(self,X):
        return set(v for v in self.descendants(X) if not self.children(v))

    def lower_boundary(self,X, context = False):
        if not context:
            closure = self.ancestors(X)
            return set(x for x in X if not (self.children(x) & closure))
        else:
            return self.lower_boundary(self.descendants(X) & context)

    def upper_boundary(self,X, context = False):
        if not context:
            closure = self.descendants(X)
            return set(x for x in X if not (self.parents(x) & closure))
        else:
            return self.upper_boundary(self.ancestors(X) & context)
    
    def lower_intersection(self, X, Y):
        return self.lower_boundary(X, X|Y) & self.lower_boundary(Y, X|Y)
    
    def upper_intersection(self, X, Y):
        return self.upper_boundary(X, X|Y) & self.upper_boundary(Y, X|Y)
    

def get_disease_phenotype_frequency_df(
    resource_dir,
    hpo_, # = RESOURCE_DIR,
    normalize_freq = True,
    fill_null = 'mean'
    # fill_null = 'median'
    # fill_null = 0
    # fill_null = float(1)
    ):

    # Read in the HPOA
    df = pd.read_csv(
        # os.path.join(resource_dir,'disease_phenotype_frequency.txt'),
        resource_dir,
        delimiter = '\t',
        skiprows = 5,
        header = None,
        low_memory = False)
    
    # Get the HPOA columns for renaming
    df.columns = ['disease_id','disease_name','qualifier','hpo_id',
              'reference','evidence','onset','frequency',
              'sex','modifier','aspect','biocuration'    
    ]

    # Limit to OMIM diseases
    df = df[df['disease_id'].str.contains('OMIM')] # Has to be an OMIM disease


    # Normalizing disease-phenotype frequencies
    if normalize_freq:
        
        df['frequency'] = df['frequency'].apply(normalize_frequency, hpo_ = hpo_)
        
        # For null values, determine how they are filled based on input.fill_null
        if type(fill_null) == float: 
            fill_freq = fill_null
            df['frequency'] = df['frequency'].fillna(fill_freq)
        elif fill_null == 'median':
            fill_freq = df['frequency'].median() # All defined frequency values
            df['frequency'] = df['frequency'].fillna(fill_freq)
        elif fill_null == 'mean': 
            fill_freq = df['frequency'].mean() # All defined frequency values
            df['frequency'] = df['frequency'].fillna(fill_freq)
        elif fill_null in [None,False]:
            pass
    
    return df  # Returns normalized HPOA


def build_propagated_frequency_map(resource_dir, hpo_): # = RESOURCE_DIR):

    # Get the HPOA with imputed/normalized frequency column (df)
    df = get_disease_phenotype_frequency_df(resource_dir, hpo_)
    
    # Intialize disease-phenotype frequency dictionary
    F = dict()

    # Loop over each row in the imputed/normalized HPOA
    for D,x,f in df[
        (df['qualifier']!='NOT') & # Not negated
        (df['aspect']=='P') # Has to be a phenotype
        ][['disease_id','hpo_id','frequency']].values:
        
        # Check to see if a key already exists, create one if it doesn't
        F[D] = F.get(D,defaultdict(float))            

        for u in hpo_.ancestors(x):
            F[D][u] = max(
                F[D].get(
                    u, # if the ancester frequency if already defined
                    0),     # if not defined, inputs 0    
                f              # if the ancestor frequency larger than current value
            )
            
    # Get rid of non-OMIM entries
    F = {k:v for k,v in F.items() if re.search('OMIM',k)}

    return F
