from Tkinter import *
import tkFileDialog
import pandas as pd
import csv,re
import scipy.stats as st
import numpy as np

root = Tk()
root.withdraw()
data_file = tkFileDialog.askopenfile(initialdir='/home/yinonbaron/Documents/Experiments')

data = np.genfromtxt(data_file.name,dtype='object', delimiter=',')
time = data[1,1:].astype('float')
data = data[3:,1:].astype('float')
norm_data = data - data[:,0:1]
window_size = 10
slopes = np.empty([shape(norm_data)[0],shape(norm_data)[1]-window_size+1])
rsquared = np.empty([shape(norm_data)[0],shape(norm_data)[1]-window_size+1])
for x,i in enumerate(norm_data):
    for j in range(0,shape(norm_data)[1]-10+1):
        res = st.linregress(time[j:j+window_size],norm_data[x,j:j+window_size])
        slopes[x,j] = res.slope
        rsquared[x,j] = res.rvalue

rates = np.max(slopes,1)