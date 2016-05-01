# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>



import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys
import csv,re
import pandas as pd

def num(s):
    if s:
        try:
            return int(s)
        except ValueError:
            return float(s)
    return -1

def findcols(s):
    list_of_cols = []
    for col in df.columns:
        print col
        if re.match(s,col):
            list_of_cols.append(col)
    return list_of_cols

def findi(s):
    list_of_indexs = []
    for i in df.index:
        label = i.split('_')[1]
        if re.match(s,label):
            list_of_indexs.append(i)
    return list_of_indexs




datafile = '/home/yinonbaron/git/prs assay/assay_251028A.csv'
mapfile = '/home/yinonbaron/git/prs assay/map_251028A.csv'



if True:
    df = pd.read_csv(datafile,index_col=0)
    dfMap = pd.read_csv(mapfile,index_col=0)

for i in df.index:
    row = i[0]
    col = i[1:]
    df.loc[i,'label'] = dfMap.loc[row,col]
df.set_index('label', inplace = True)


# <codecell>

#df.loc[findi('n1b')].transpose().plot(figsize=(10,10))

df.loc[findi('hxa')]