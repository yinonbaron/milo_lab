from Tkinter import *
import tkFileDialog
import pandas as pd
import csv,re
import scipy.stats as st

def num(s):
    if s:
        try:
            return int(s)
        except ValueError:
            return float(s)
    return -1

def findcols(s):
    list_of_cols = []
    for col in data.columns:
        print col
        if re.match(s,col):
            list_of_cols.append(col)
    return list_of_cols

def findi(s,df):
    list_of_indexs = []
    for i in df.index:
        '''        
        label = i.split('_')[1]
        if re.match(s,label):
            '''
        if re.findall(s,i):
            list_of_indexs.append(i)
    return list_of_indexs


root = Tk()
root.withdraw()
data_file = tkFileDialog.askopenfile(initialdir='/home/yinonbaron/Documents/Experiments')
map_file = tkFileDialog.askopenfile(initialdir='/home/yinonbaron/Documents/Experiments')

data = pd.DataFrame.from_csv(data_file,index_col=0)
map_df = pd.DataFrame.from_csv(map_file,index_col=0)
time = data.loc[data.index[0]]
data = data.loc[data.index[2:]]

for i in data.index:
    row = i[0]
    col = i[1:]
    data.loc[i,'label'] = map_df.loc[int(col),row]
data.set_index('label', inplace = True)

norm_data = data.subtract(data.icol(0),axis=0)

window_size = 10
slopes = np.empty([len(norm_data.index),len(data.columns)-10+1])
rsquared = np.empty([len(norm_data.index),len(data.columns)-10+1])
for x,i in enumerate(norm_data.index):
    for j in range(0,len(data.columns)-10+1):
        res = st.linregress(time[j:j+window_size],norm_data.loc[i][j:j+window_size])
        slopes[x,j] = res.slope
        rsquared[x,j] = res.rvalue

slopes_df = pd.DataFrame(slopes,index = norm_data.index)
rvalue_df = pd.DataFrame(rsquared,index = norm_data.index)

rates = pd.DataFrame.max(slopes_df,axis=1)
