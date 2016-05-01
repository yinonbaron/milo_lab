# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:33:53 2015

@author: yinonbaron
"""
from models import *
from collections import Counter
import uncertainties.unumpy as unumpy  
from model_addons import add_to_model
from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
from collections import defaultdict
from thermodynamics_for_cobra import reaction_thermodynamics
model = create_cobra_model_from_sbml_file('/home/yinonbaron/git/shared_data/iJO1366.xml')
convert_to_irreversible(model)
add_to_model(model)
tr = reaction_thermodynamics(model.reactions)
cofactors = ['10fthf_c','mlthf_c','h_c','adp_c','atp_c','nad_c','h2o_c','nadh_c','pi_c','nadp_c','nadph_c','coa_c','accoa_c','pyr_c','h2o2_c','ppi_c','gtp_c','h_p','h2o_p','q8h2_c','mql8_c','q8_c','o2_c','gmp_c','amp_c','h2s_c','mqn8_c','glu_L_c','gln_L_c','co2_c','nh4_c','pi_p','ppgpp_c','gdp_c','fmn_c','so3_c','gthox_c','zn2_p']
del_g_tresh = -20
metabolite_set = list(model.metabolites)

    
def create_click(metabolite):
    click = list()
    #click.append(metabolite)
    metabolite_set.pop(metabolite_set.index(metabolite))
    next_step = [metabolite]
    while len(next_step) > 0:
        for item in next_step:
            click.append(item)
        next_step = advance_step(next_step)
        '''
    neighbors = reversible_neigbors(metabolite)
    if len(neighbors) >0:
        click.append(neighbors)
        next_step = neighbors
        next_neighbors = list()
        for item in next_step:
            next_neighbors.append(reversible_neigbors(metabolite))
        next_neighbors = np.ndarray.tolist(np.reshape(np.array(next_neighbors),[1,-1]))[0]
        while len(neighbors) >0:
            neighbors = reversible_neigbors(metabolite)
            if len(neighbors) >0:
                click.append(neighbors)
            print click
            '''
    return click

def advance_step(metabolite_list):
    next_neighbors = list()
    for item in metabolite_list:
        s = reversible_neigbors(item)
        if len(s)>0:
            next_neighbors.append(s)
    next_neighbors = [item for sub in next_neighbors for item in sub]
    next_neighbors = list(set(next_neighbors)) # remove dups
    if len(next_neighbors) > 0:
        for item in next_neighbors:
            metabolite_set.pop(metabolite_set.index(item))
    return next_neighbors

def reversible_neigbors(metabolite):
    neighbor_reactions = [r for r in model.reactions if metabolite in r.reactants]
    reversible_neighbor_reactions = [r for r in neighbor_reactions if r.dGm_prime > del_g_tresh]
    reversible_neighbor_metabolites = [r.products for r in reversible_neighbor_reactions]
    reversible_neighbor_metabolites = [item for sub in reversible_neighbor_metabolites for item in sub]
    reversible_neighbor_metabolites = [r for r in reversible_neighbor_metabolites if r in metabolite_set]
    return reversible_neighbor_metabolites
    

clicks = list()
for cof in cofactors:
    ind = [r for r in metabolite_set if cof == r.id]
    #print ind
    if len(ind)>0:
        metabolite_set.pop(metabolite_set.index(ind[0]))

while len(metabolite_set) > 0:
    metabolite = metabolite_set[0]
    clicks.append(create_click(metabolite))
    
x = list()
for i in clicks:
    x.append(len(i))