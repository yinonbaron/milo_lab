#!/usr/bin/python
import scipy.sparse
from itertools import chain, combinations
import os, sys, pickle, csv
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

from analysis_toolbox import model_summary, plot_multi_PPP
from models import *
from optknock import OptKnock
from draw_flux import DrawFlux
from html_writer import HtmlWriter

NUMBER_OF_KNOCKOUTS = 1
BIOMASS_PRECURSORS = set(['3pg', 'accoa', 'e4p', 'f6p', 'g3p', 'g6p', 'gln_L', 'glu_L',
                      'oaa', 'pep', 'pyr', 'r5p'])

single_ko_list = ['ICL', 'GLCpts','PGI','PFK','FBP','FBA','TPI','GAPD,PGK','PYK','PPS','PDH','PFL',
                 'GND','EDD,EDA','RPE','RPI','TKT1,TKT2','TALA','PPC',
                 'PPCK','CS','ACONTa,ACONTb','ALCD2x','ACALD']

single_ko_list = ['PFK,PGM,TALA,ICL']

carbon_sources = ['g6p', 'f6p', '6pgc', 'r5p', 'succ', 'xu5p_D', '2pg', 'pyr', 'ac', 'dhap']

carbon_sources = ['xu5p_D']

def generate_model(ko_list, carbon_source,
                   carbon_uptake_rate=50, BM_lower_bound=0.1,
                   knockins='EDD,EDA,RBC,PRK,SBP,SBA'):
    """
        carbon_uptake_rate is in [mmol C / (gDW*h)]
    """
    model = init_wt_model('core', {}, BM_lower_bound=BM_lower_bound)
    knockin_reactions(model, knockins, 0, 1000)
    knockin_reactions(model, 'EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,EX_2pg,EX_e4p,EX_6pgc', 0, 0)
    knockout_reactions(model, 'G6PDH2r,PGL') # always knockout the ZWF gene
    
    if carbon_source == 'electrons':
        knockin_reactions(model, 'RED', 0, carbon_uptake_rate*2)
    else:
        # find out how many carbon atoms are in the carbon source
        # and normalize the uptake rate to be in units of mmol carbon-source / (gDW*h) 

        nC = 0
        for cs in carbon_source.split(','):
            met = model.metabolites[model.metabolites.index(cs + '_c')]
            nC += met.formula.elements['C']
        uptake_rate = carbon_uptake_rate / float(nC)

        for cs in carbon_source.split(','):
            set_exchange_bounds(model, cs, lower_bound=-uptake_rate)
    
    for ko in ko_list:
        knockout_reactions(model, ko)
    return model

def get_precursors(ko_list, energy_source='electrons', knockins='EDD,EDA,RBC,PRK'):
    """
        checks whether the model can produce at least one biomass precursor
    """
    temp_model = generate_model(ko_list, energy_source, BM_lower_bound=0, knockins=knockins)
    electron_precursors = set([])
    for precursor in BIOMASS_PRECURSORS:
        set_single_precursor_objective(temp_model, precursor)

        # Run the optimization for the objective reaction and medium composition
        # set in the file.
        ok = OptKnock(temp_model)
        ok.prepare_FBA_primal()
        ok.solve()

        if ok.get_objective_value() > 1e-5:
            electron_precursors.add(precursor)
    return electron_precursors     

def analyze_rubisco_dependent():
    
    target_reaction = 'RBC'
    
    print "There are %d single knockouts\n" % len(single_ko_list)
    print "There are %d carbon sources: %s\n" % (len(carbon_sources), ', '.join(carbon_sources))
    
    csv_out = csv.writer(open('res/semi_autotrophic.csv', 'w'))
    # csv header    
    csv_out.writerow(['trophism', 
                      'ko_list', 
                      'unique_carbon_precursors', 
                      'unique_electrons_precursors', 
                      'carbon_source'])
                      
    for i in xrange(1, NUMBER_OF_KNOCKOUTS+1):
        for ko_list in combinations(single_ko_list, i):            
            print ko_list

            electron_precursors = get_precursors(ko_list, energy_source='electrons')
            # check if KO strain has the potential to be an autotroph
            if electron_precursors:
                trophism = 'potential'
                print 'semi-autotrophic potential'
                
                for carbon_source in carbon_sources:
                    temp_model = generate_model(ko_list, carbon_source)
                    biomass_yield = OptKnock(temp_model).solve_FBA() or np.nan
                    if np.isnan(biomass_yield):
                        # find carbon source biomass precursors
                        continue                        
                    slope = OptKnock(temp_model).get_slope(target_reaction)
                    if slope > 0:
                        print slope
                        # if all biomass precursors can be generated using electrons this in an autotroph    
                        if electron_precursors == BIOMASS_PRECURSORS:
                            trophism = 'autotroph'
                            print "autotroph - not yuck!"
                            carbon_precursors = []
                            csv_out.writerow([trophism, 
                                            ';'.join(sorted(ko_list)), 
                                            ';'.join(carbon_precursors), 
                                            ';'.join(electron_precursors)])
                            break
                        else: # at least one elctron precursor is synthesized via "carbon-fixation", thus may be a semi
                            # remove RBC from model and grow on carbon to find BM precursors generated by carbon only
                            temp_model = generate_model(ko_list, carbon_source, BM_lower_bound=0, knockins='EDD,EDA')
                            carbon_precursors = set([])
                            for precursor in BIOMASS_PRECURSORS:
                                set_single_precursor_objective(temp_model, precursor)
                                # Run the optimization for the objective reaction and medium composition set in the file.
                                ok = OptKnock(temp_model)
                                ok.prepare_FBA_primal()
                                ok.solve()
                                if ok.get_objective_value() > 1e-5:    
                                    carbon_precursors.add(precursor)
                                    print precursor
    
                            if set(BIOMASS_PRECURSORS) - carbon_precursors <= electron_precursors:
                                trophism = 'semi-autotroph'
                                print "this is a semi-autotroph! can grow noramlly on:\t", carbon_source
                            csv_out.writerow([trophism, 
                                            ';'.join(sorted(ko_list)), 
                                            ';'.join(carbon_precursors), 
                                            ';'.join(electron_precursors), 
                                            carbon_source])

                            

def main():
    analyze_rubisco_dependent()
    #print test_semi_potential('', '')
    return
    
#    temp_model = generate_model('', '', 'g6p', BM_lower_bound=0.1)
#    ok = OptKnock(temp_model)
#    print ok.solve_FBA()
#    print ok.solution.x
#    main_html = HtmlWriter('res/semi_auto.html')
#    ok.draw_svg(main_html)
#    ok.model_summary(main_html)
#    return
#
#    main_html.write('<h1>Flux Balance Analysis of semi-autotrophic E. coli</h1>\n')
#    
#    main_html.write('<table border="1">\n')
#    main_html.write('<ul>\n')
#    for precursor in BIOMASS_PRECURSORS:
#        model = init_wt_model('core', {}, BM_lower_bound=0)
#        knockout_reactions(model, 'PGM')
#        knockin_reactions(model, 'RED,RBC,PRK')
#        
#        set_single_precursor_objective(model, precursor)
#
#        # Run the optimization for the objective reaction and medium composition
#        # set in the file.
#        ok = OptKnock(model)
#        ok.prepare_FBA_primal()
#        ok.solve()
#        growth = ok.get_objective_value()
#        print '%s: f = %.3g' % (precursor, growth)
#        main_html.write('<li>%s Yield: %.3g</li>\n' % (precursor, growth))
#    main_html.write('</ul>\n')
#
#        #if growth is None:
#        #    main_html.write('<h2>No growth possible</h2>\n')
#        #else:
#        #    print ': f = %.3g' % growth
#        #    main_html.write('<h2>Growth Yield: %.3g</h2>\n' % growth)
#        #    ok.draw_svg(main_html)
#        #    ok.model_summary(main_html)
#    main_html.close()

if __name__ == "__main__":
    main()
