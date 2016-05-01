# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 17:51:39 2015

@author: yinonbaron
"""

import scipy.sparse
#from itertools import chain, combinations
import os, sys, pickle, csv
#from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

#from analysis_toolbox import model_summary, plot_multi_PPP
from models import *
#from optknock import OptKnock
#from draw_flux import DrawFlux
#from html_writer import HtmlWriter

#BM_lower_bound=0.1
#model = init_wt_model('full', {}, BM_lower_bound=BM_lower_bound)
#add_metabolite_exchange(model, 'r5p', -50, upper_bound=0)
#model.reactions[model.reactions.index('EX_r5p_e')].lower_bound = -50
#model.reactions[model.reactions.index('EX_r5p_e')].upper_bound = 0
#set_single_precursor_objective(model,'atp')
#solution = model.optimize()
#flux = solution.x
#t = model.reactions[model.reactions.index('Biomass_atp')]
#print t.objective_coefficient

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
r5p = model.metabolites.get_by_id('r5p_c')
prpp = model.metabolites.get_by_id('prpp_c')
water = model.metabolites.get_by_id('h2o_c')
r5p_consuming = [r for r in model.reactions if r5p in r.reactants]
pyruvate = model.metabolites.get_by_id('pyr_c')
oaa = model.metabolites.get_by_id('oaa_c')
asp =model.metabolites.get_by_id('asp_L_c')
aca =model.metabolites.get_by_id('accoa_c')


prpp_consuming = [r for r in model.reactions if prpp in r.reactants]
pur_consuming = [r for r in model.reactions if pyruvate in r.reactants]
oaa_consuming = [r for r in model.reactions if oaa in r.reactants]
asp_consuming = [r for r in model.reactions if asp in r.reactants]
aca_consuming = [r for r in model.reactions if aca in r.reactants]

irreversible_reactions = [r for r in model.reactions if r.dGm_prime < -10]
result = dict()#defaultdict(list)
for x,reaction in enumerate(irreversible_reactions):
    tmp = dict()
    for metabolite in reaction.products:
        if metabolite.id not in ['h_c','adp_c','atp_c','nad_c','h2o_c','nadh_c','pi_c','nadp_c','nadph_c','coa_c','accoa_c','pyr_c','h2o2_c','ppi_c','gtp_c','h_p','h2o_p','q8h2_c','mql8_c','q8_c','o2_c','gmp_c','amp_c','h2s_c','mqn8_c','glu_L_c','gln_L_c','co2_c','nh4_c','pi_p','ppgpp_c','gdp_c','fmn_c','so3_c','gthox_c']:
            next_node = [r for r in model.reactions if metabolite in r.reactants]
            tmp[metabolite] = [r for r in next_node if r.dGm_prime < -10] #tmp[metabolite] = [r for r in next_node if r.logRI > 5]
            result[reaction] = tmp # result[x].append({metabolite: [r for r in next_node if r.dGm_prime < -30]})

sumres = [len([item for sub in w.values() for item in sub]) for w in result.itervalues()]
q = np.squeeze(np.array([[[result.keys()[x],result[result.keys()[x]].keys()]] for x,r in enumerate(sumres) if r>0]))
filt_result = {key: result[key] for key in [result.keys()[x] for x,r in enumerate(sumres) if r>0]}
p = [len([r for r in model.reactions if filt_result[filt_result.keys()[x]].keys()[0] in r.reactants]) for x,s in enumerate(filt_result.keys())]
'''
sumres = list()
for i in range(0,len(result)):
    record = result[i]
    sumres.append(sum([len(metabolite.values()[0]) for metabolite in record]))
    '''
[x for x,r in enumerate(sumres) if r>7]
data = np.empty([len(np.where(np.array(p)==1)[0]),6],dtype='object')
enum = 0;
for x,record in enumerate(p):
    if record == 1:
        for meta in filt_result[filt_result.keys()[x]].iteritems():
            if len(meta[1])>0:
                data[enum,:] = [filt_result.keys()[x].subsystem,filt_result.keys()[x].name, filt_result.keys()[x].gene_name_reaction_rule,meta[0].name,[r.name for r in meta[1]][0],[r.gene_name_reaction_rule for r in meta[1]][0]]
        enum = enum+1
np.savetxt('IRPairs.tsv',data,'%s','\t')

            
'''
reaction_log = list()
reaction_log.append(prpp_consuming)
for n in range(1,3):
    for r in reaction_log[n-1]:
        tmp = list()
        tmp.append({r: [q for q in model.reactions if r in q.reactants]})
    reaction_log.append(tmp)        
''' 


