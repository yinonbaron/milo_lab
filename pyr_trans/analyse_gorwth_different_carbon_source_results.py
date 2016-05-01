# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
#import  seaborn as sns

import matplotlib.pyplot as plt

#current_palette = sns.color_palette()
#sns.palplot(current_palette)
#sns.palplot(sns.color_palette("hls", 8))


e = pd.read_csv('/home/yinonbaron/Documents/Experiments/pyruvate transporter/20151210_F-pyr_resistance_colonies_different_cs_partial.csv')
time = e['Time [s]']
data = e[e.columns[5:]]


e2 = pd.read_csv('/home/yinonbaron/Documents/Experiments/pyruvate transporter/20151220_pyr_trans_clones_different_cs.csv')
time2 = e2.loc[0][e2.columns[1:]]
data2 = e2.loc[2:][e2.columns[1:]]


raw_data = np.loadtxt('/home/yinonbaron/Documents/Experiments/pyruvate transporter/20151220_pyr_trans_clones_different_cs.csv',delimiter=',',dtype='object')

time = raw_data[1,1:].astype('float')
data = raw_data[3:,1:].astype('float')
reshaped_data = np.reshape(data,[8,12,-1])
norm_data = np.subtract(reshaped_data, reshaped_data[:,:,0:1])

plt.plot(time[:,], norm_data[0,0:1,:].T)