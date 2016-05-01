import json
import pandas as pd
import xmltodict
import numpy as np
def get_record_info(record,query_len):
    res = np.empty([4,],dtype='object')
    if type(record['Hit_hsps'].values()[0]) == list:
            start = min([float(r['Hsp_query-from']) for r in record['Hit_hsps'].values()[0]])
            end = max([float(r['Hsp_query-to']) for r in record['Hit_hsps'].values()[0]])
            res[0] = record['Hit_accession']
            res[1] = (end-start+1)/float(query_len)
            res[2] = max([float(r['Hsp_identity'])/(float(r['Hsp_query-to'])-float(r['Hsp_query-from'])+1) for r in record['Hit_hsps'].values()[0]])
            res[3] = min([float(r['Hsp_evalue']) for r in record['Hit_hsps'].values()[0]])           
    else:
            # Just one hsps
            start = float(record['Hit_hsps']['Hsp']['Hsp_query-from'])
            end = float(record['Hit_hsps']['Hsp']['Hsp_query-to'])
            res[0] = record['Hit_accession']
            res[1] = (end-start+1)/float(query_len)
            res[2] = float(record['Hit_hsps']['Hsp']['Hsp_identity'])/(end-start+1)
            res[3] = float(record['Hit_hsps']['Hsp']['Hsp_evalue'])
    return(res)

tmhmm = pd.read_csv('tmhmm_results.csv',index_col=0)

json_file = open('/home/yinonbaron/git/pubique_tm_genes/6GH5TK4Z015-Alignment.json').read()
data = json.loads(json_file)

with open('/home/yinonbaron/Downloads/6JF3W81A015-Alignment.xml') as fd:
    obj = xmltodict.parse(fd.read())

gene_names = [r['Iteration_query-def'][:r['Iteration_query-def'].find(' ')] for r in obj['BlastOutput']['BlastOutput_iterations']['Iteration']]
align_df = pd.DataFrame(index= gene_names,columns = ['Match','Coverage','Identity','E-val'])
for i in range(0,len(gene_names)):
    record = obj['BlastOutput']['BlastOutput_iterations']['Iteration'][i]
    if record['Iteration_hits'] == None:
        continue;
    if type(record['Iteration_hits'].values()[0]) == list:
#        continue;
        tmp = np.empty([4,len(record['Iteration_hits'].values()[0])],dtype='object')            
        for x,hit in enumerate(record['Iteration_hits'].values()[0]):
            tmp[:,x] = get_record_info(hit,record['Iteration_query-len']).T
            best_hit = tmp[1,:] == max(tmp[1,:])
            align_df['Match'][gene_names[i]] = tmp[0,best_hit][0]
            align_df['Coverage'][gene_names[i]] = tmp[1,best_hit][0]
            align_df['Identity'][gene_names[i]] = tmp[2,best_hit][0]
            align_df['E-val'][gene_names[i]] = tmp[3,best_hit][0]
    else:
        tmp = get_record_info(record['Iteration_hits']['Hit'],record['Iteration_query-len'])
        align_df['Match'][gene_names[i]] = tmp[0]
        align_df['Coverage'][gene_names[i]] = tmp[1]
        align_df['Identity'][gene_names[i]] = tmp[2]
        align_df['E-val'][gene_names[i]] = tmp[3]
        # Just one hit
#        if type(record['Iteration_hits']['Hit']['Hit_hsps'].values()[0]) == list:
#            start = min([float(r['Hsp_query-from']) for r in record['Iteration_hits']['Hit']['Hit_hsps'].values()[0]])
#            end = max([float(r['Hsp_query-to']) for r in record['Iteration_hits']['Hit']['Hit_hsps'].values()[0]])
#            align_df['Match'][gene_names[i]] = record['Iteration_hits']['Hit']['Hit_accession']
#            align_df['Coverage'][gene_names[i]] = (end-start+1)/float(record['Iteration_query-len'])
#            align_df['Identity'][gene_names[i]] = max([float(r['Hsp_identity'])/(float(r['Hsp_query-to'])-float(r['Hsp_query-from'])+1) for r in record['Iteration_hits']['Hit']['Hit_hsps'].values()[0]])
#            align_df['E-val'][gene_names[i]] = min([float(r['Hsp_evalue']) for r in record['Iteration_hits']['Hit']['Hit_hsps'].values()[0]])           
#        else:
            # Just one hsps
#            start = float(record['Iteration_hits']['Hit']['Hit_hsps']['Hsp']['Hsp_query-from'])
#            end = float(record['Iteration_hits']['Hit']['Hit_hsps']['Hsp']['Hsp_query-to'])
#            align_df['Match'][gene_names[i]] = record['Iteration_hits']['Hit']['Hit_accession']
#            align_df['Coverage'][gene_names[i]] = (end-start+1)/float(record['Iteration_query-len'])
#            align_df['Identity'][gene_names[i]] = float(record['Iteration_hits']['Hit']['Hit_hsps']['Hsp']['Hsp_identity'])/(end-start+1)
#            align_df['E-val'][gene_names[i]] = float(record['Iteration_hits']['Hit']['Hit_hsps']['Hsp']['Hsp_evalue'])

uni_df = pd.read_csv('uniprot_mapping.csv',sep='\t',index_col=2)
merged_df = align_df.merge(uni_df, left_on= align_df['Match'],right_index=True)
tm_merged_df = merged_df.join(tmhmm,how='inner')

