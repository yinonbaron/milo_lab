# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:21:35 2016

@author: yinonbaron
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

od = pd.read_csv('/home/yinonbaron/Downloads/od.csv',sep=',',)
flour = pd.read_csv('/home/yinonbaron/Downloads/flour.csv',sep=',')
time_list = pd.read_csv('/home/yinonbaron/Downloads/time.csv',sep=',')
time_list = time_list.columns.astype('float')
flour.replace('OVER',np.nan,inplace=True)

od_np = od.values.reshape([8,12,-1])
f_np = flour.values.reshape([8,12,-1]).astype('float')

bg = od_np[:,0,:].mean(axis=1).reshape([1,-1])
f_bg = f_np[:,0,:].mean(axis=1).reshape([1,-1])

od_np_sub_bg = (np.swapaxes(od_np,0,2)-bg).swapaxes(0,2)
f_np_sub_bg = (np.swapaxes(f_np,0,2)-f_bg).swapaxes(0,2)


f_norm_od = f_np_sub_bg/od_np_sub_bg

sig1 = (f_norm_od[1:6,10,:]-f_norm_od[1:6,11,:])/(f_norm_od[1:6,8,:]-f_norm_od[1:6,11,:])
sig2 = (f_norm_od[1:6,6,:]-f_norm_od[1:6,2,:])/(f_norm_od[1:6,4,:]-f_norm_od[1:6,2,:])
activity1 = list()
activity2 = list()
for i in range(1,6):
    activity1.append(np.nanmean(sig1[i-1,(od_np_sub_bg[i,10,:] > 0.1) & (od_np_sub_bg[i,10,:] < 0.2)]))
    activity2.append(np.nanmean(sig2[i-1,(od_np_sub_bg[i,6,:] > 0.1) & (od_np_sub_bg[i,6,:] < 0.2)]))
e = list()
od_pd_sub_bg = pd.DataFrame(od_np_sub_bg.reshape([96,-1]).T)
for i in od_pd_sub_bg.columns:
    model = pd.ols(y = np.log(od_pd_sub_bg[i]),x= pd.Series(time_list/3600),window_type='rolling',window=5)
    e.append(model.beta['x'])

gr = pd.DataFrame(e).values.reshape([8,12,-1])

plt.plot(gr[1:6,10,:].max(axis=1),activity1,'*')
plt.figure()
plt.plot(gr[1:6,6,:].max(axis=1),activity2,'*')