import pandas as pd
candidates = ['yadH','ycjP','yegH','yciB','yigM','ychE','yqaA','ycaD','yeiH','yhdX','yhdY','ydiN','ybhL','yhhQ']
gene_list = pd.read_csv('ecoli_genes.csv',index_col=1)


expression_data = pd.read_csv('/home/yinonbaron/git/proteomics-collection/copies_fL/ecoli_Schmidt_et_al_2015.csv',index_col=0)
x = expression_data.loc[gene_list.loc[candidates]['Locus Name']]

a = expression_data['GLC_BATCH_mu=0.58_S']
a.replace([np.inf,-np.inf,np.nan,0],0.00001,inplace=True)
b = (expression_data['PYR_BATCH_mu=0.40_S']/a)
b2 = b.loc[b>0]

e = gene_list.merge(pd.DataFrame(b2.loc[b2>b2.quantile(.95)]),how='inner',left_on='Locus Name',right_index=True)
e.sort(columns=0,inplace=True)
f = e.merge(expression_data[['PYR_BATCH_mu=0.40_S','GLC_BATCH_mu=0.58_S']],how='inner',left_on='Locus Name',right_index=True)
