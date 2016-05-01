# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 17:25:45 2015

@author: yinonbaron
"""

## Labels
def add_labels(labels,label_size,start_location,direction,label_mat):
    for j in start_location:
        for i,x in enumerate(labels):
            label_tmp = tile(x,label_size)
            if direction == 'down':
                label_mat[j[0]+i:j[0]+label_size[0]*i,j[1]:j[1]+label_size[1]] = label_tmp
            elif direction == 'side':
                label_mat[j[0]:j[0]+label_size[0],j[1]+i:j[1]+label_size[1]*i] = label_tmp
    return label_mat

'''
        if direction == 1:
            if x*start_location[1]+<=12:
                
                label[i] = [[start_location[0],start_location[1]+x*label_size[1]],[start_location[0]+,start_location[1]+(x+1)*label_size[1]]]
            else:
                
        if direction == 2:
            label[i] = [[start_location[0]+x*label_size[0],start_location[1]+x*label_size[1]],[start_location[0]+(x+1)*label_size[0],start_location[1]+(x+1)*label_size[1]]]
            
def insert_side():
    
def insert_down():
    

#def generate_map(labels,outfile):
import numpy as np
labels = {};
concentrations = [0,0.05,0.075,0.1,0.15,0.2,0.4,1]
for i, x in enumerate(concentrations):
    label[str(i)] = [[x+1,1],[x+1,12]]
concentration_labels = tile(np.array(concentrations,dtype='str'),[12,1]).transpose()
concentration_locs = np.empty([12,8])
#concentration_locs[:,:] = concentration_labels


#strain_labels

strains = ['wt','n1','n2a','n3']
'''


concentration= [0,0.05,0.75,0.1,0.15,0.2]
conc_label = add_labels(concentration,label_size = [8,1],start_location = [[0,0],[0,6]],direction='side', label_mat =np.empty([8,12]))

strains = ['wt','n1','n2','n3']
strain_label = add_labels(strains,label_size = [4,6],start_location = [[0,0]],direction='side', label_mat =np.empty([8,12],dtype='object'))