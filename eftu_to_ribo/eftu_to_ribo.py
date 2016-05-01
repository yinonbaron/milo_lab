# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 14:27:44 2016

@author: yinonbaron
"""

import pandas as pd
import matplotlib.pyplot as plt
copy_fL = pd.read_csv('../copies_fL/ecoli_Schmidt_et_al_2015.csv',sep=',',index_col=0)

eftu = ['b3339','b3980']
ribo = ['b2330','b3984','b3317','b3320','b3319','b3308','b3305','b4203','b3985','b3983','b3986','b3231','b3310','b3301','b3313','b3294','b3304','b2606','b1716','b3186','b3315','b3318','b3309','b2185','b3185','b3637','b3312','b3302','b3936','b1089','b3636','b3703','b1717','b3299','b0911','b0169','b3314','b3296','b3303','b4200','b3341','b3306','b3230','b3321','b3297','b3342','b3298','b3307','b3165','b2609','b3311','b4202','b3316','b0023','b3065','b1480']
eftu_ab = copy_fL.ix[eftu].sum()
ribo_ab = copy_fL.ix[ribo].sum()/len(ribo)

gr = [float(r.split('=')[1].split('_')[0]) for r in ribo_ab.index.tolist()]
ax = plt.figure()
plt.plot(gr,eftu_ab/ribo_ab,'.',figure = ax)
plt.ylim(0,20)
plt.xlim(0,2.5)