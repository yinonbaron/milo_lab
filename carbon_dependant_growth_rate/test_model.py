# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 15:57:46 2015

@author: yinonbaron
"""

from cobra.io.sbml import create_cobra_model_from_sbml_file

#model = create_cobra_model_from_sbml_file("../shared_data/ecoli_core_model.xml")
model = create_cobra_model_from_sbml_file("../cobrapy/cobra/test/data/iJO1366.xml")
# the core model has these annoying '_b' metabolites that are used as
# 'ghost' metabolites that balance the exchange reactions. they should
# be ignored in the mass-balance equation and therefore the best way to
# deal with them is to remove them from all the reactions
for m in model.metabolites:
    if m.id.endswith('_b'):
        for r in m.reactions:
            coeff = r.get_coefficient(m)
            r.add_metabolites({m : -coeff})
rxns = {r.id:r for r in model.reactions}
rxns['Biomass_Ecoli_core_w_GAM'].upper_bound = 0.89
rxns['EX_glc_e'].lower_bound = -14
rxns['EX_ac_e'].upper_bound = 7.5
rxns['EX_o2_e'].lower_bound = -19.8
solution = model.optimize()
print solution.f