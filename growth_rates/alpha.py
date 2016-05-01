# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 10:24:59 2016

@author: yinonbaron
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize.minpack import curve_fit

file_name = '/home/yinonbaron/Downloads/20160331_pyr_fpyr_mic_full.csv'

x = pd.read_csv(file_name)
#x[x.ix[:,0].str.find('Time [s]')==0].
x2 = x.ix[1:,1:]
x2 = x2.transpose()
x3 = x2.reset_index()
x3['index'] = x3['index'].astype('float')
blank_wells = [86,87,88,89]
x4 = x3
x4.ix[:,1:] = x4.ix[:,1:] -x3.iloc[:,blank_wells].as_matrix().mean()
x4.ix[:,1:] = np.log2(x4.ix[:,1:])
x5 = x4.iloc[:,~pd.isnull(x4.sum()).as_matrix()]
#time = x.ix[1:,1:].columns.astype('float')
e = list()
for i in x5.columns[1:]:
    model = pd.ols(y = x5[i],x=x5['index']/3600,window_type='rolling',window=30)
    e.append(model.beta['x'])
    
t = pd.DataFrame(e,index=x5.columns[1:])
mean_slopes = t.iloc[:,:40].transpose().mean()
mean_slopes = mean_slopes.iloc[:-8].as_matrix()
result_slopes = mean_slopes.reshape([-1,10])
conc = [0,0.2/128,0.2/64,0.2/32,0.2/16,0.2/8,0.2/4,0.2/2,0.2]
conc.reverse()
plt.errorbar(conc,result_slopes.mean(axis=0)[:-1],result_slopes.std(axis=0)[:-1])

guess_a, guess_b, guess_c = 0.2, -6, 0.01
guess = [guess_a, guess_b, guess_c]

exp_decay = lambda x, A, t, y0: A * np.exp(x * t) + y0

params, cov = curve_fit(exp_decay, conc, result_slopes.mean(axis=0), p0=guess)

A, t, y0 = params

best_fit = lambda x: A * np.exp(t * x) + y0
plt.plot(conc, best_fit(np.array(conc)), 'r')
plt.xlabel('Fluoropyruvate concentration (mM)')
plt.ylabel('Growth rate (h^-1)')