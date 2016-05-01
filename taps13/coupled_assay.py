# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 13:49:01 2015

@author: dan
"""

from collections import Counter
import uncertainties.unumpy as unumpy  
from model_addons import add_to_model
from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
from collections import defaultdict
from thermodynamics_for_cobra import reaction_thermodynamics

def filter_reactions(metabolite, only_irreversible=True):
    metabolite_forming_reactions = [r for r in metabolite.reactions 
                              if metabolite in r.products]
    
    homomeric_reactions = filter(lambda r: 'and' not in r.gene_reaction_rule, 
                                 metabolite_forming_reactions)
    reactions = {r:r.dG0_prime for r in homomeric_reactions}
    if only_irreversible:
        reactions = {k:v for k,v in reactions.iteritems() if v<10}
    return reactions

#model = create_cobra_model_from_sbml_file('../data/iJO1366.xml')
model = create_cobra_model_from_sbml_file('/home/yinonbaron/git/shared_data/iJO1366.xml')

add_to_model(model)
tr = reaction_thermodynamics(model.reactions)
nadh = model.metabolites.get_by_id('nadh_c')    
# reactions which produce NADH, and are thermodynamically favorable
# Also, all reactions are catalyzed by homomeric enzymes.
nadh_forming = filter_reactions(nadh, only_irreversible=True)

water = model.metabolites.get_by_id('h2o_c')

first_hop = defaultdict(list)
for r in nadh_forming:
    for m in r.reactants:
        if m not in [water, nadh]:
            homomeric_and_favorable = filter_reactions(m)
            for k in homomeric_and_favorable.iterkeys():
                if k not in nadh_forming:
                    first_hop[k].append([r.id,m.id])

second_hop = defaultdict(list)
for r in first_hop.keys():
    for m in r.reactants:
        if m not in [water, nadh]:
            homomeric_and_favorable = filter_reactions(m)
            for k in homomeric_and_favorable.iterkeys():
                if k not in first_hop and k not in nadh_forming:
                    second_hop[k].append([r.id,m.id])

enzyme_assays = nadh_forming.keys()+first_hop.keys()+second_hop.keys()
reaction_counter = Counter([r.subsystem for r in model.reactions])
tmp = Counter([r.subsystem for r in enzyme_assays])
coverage = {k:tmp[k]/float(reaction_counter[k])*100 for k in reaction_counter.iterkeys()}
print sum(coverage.values()) / len(coverage.keys())

#homomeric_reactions = defaultdict(list)
#reactions = []
#for r in nadh_forming_reactions:
#    g = list(r.genes)[0]
#    homomeric_reactions[r.id].append(r.name)
#    homomeric_reactions[r.id].append(g.id)
#    homomeric_reactions[r.id].append(g.name)
#    homomeric_reactions[r.id].append(g.MW)
#    homomeric_reactions[r.id].append(r.kegg_reaction.write_formula())
#    homomeric_reactions[r.id].append(unumpy.nominal_values(r.dG0_prime))
#    homomeric_reactions[r.id].append(unumpy.nominal_values(r.dGm_prime))
#    homomeric_reactions[r.id].append(unumpy.std_devs(r.dGm_prime))
##homomeric_reactions[r.id].append
#out = pd.DataFrame.from_dict(homomeric_reactions, 'index')
#out.columns = ['reaction_name','bnumber','gene_name','MW[Da]','reaction','dG0','dGm','dG_std']
#out.to_csv('../res/nadh_forming_enzymes.tsv', sep='\t')