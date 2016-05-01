# -*- coding: utf-8 -*-

import keggrest

def reduce_list(list_of_lists):
    list = [item for sublist in list_of_lists for item in sublist]
    return list

r5p = keggrest.keggrest.KEGGget('C00117')

reactions = reduce_list([r.split(' ') for r in r5p['REACTION']])
for reaction in reactions:
    keggrest.keggrest.KEGGget('C00117')
