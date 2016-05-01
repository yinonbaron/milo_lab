# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 14:25:45 2016

@author: yinonbaron
"""

from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible, revert_to_reversible
#model = create_cobra_model_from_sbml_file("../cobrapy/cobra/test/data/iJO1366.xml")
model = create_cobra_model_from_sbml_file("../shared_data/ecoli_core_model.xml")
model.reactions.Biomass_Ecoli_core_w_GAM.delete()
convert_to_irreversible(model)
num_metabolites = len(model.metabolites)
num_reactions = len(model.reactions)
total_list = model.reactions + model.metabolites
nodes_mapping = dict()
map(lambda k, v: nodes_mapping.update({k: v}), total_list, range(0,num_metabolites+num_reactions))

network = list()

for reaction in model.reactions:
    for reactant in reaction.reactants:
        network.append((nodes_mapping[reactant],nodes_mapping[reaction]))
    for product in reaction.products:
        network.append((nodes_mapping[reaction],nodes_mapping[product]))
        
n= num_metabolites + num_reactions
E= network

blocked = set([])
B={}
AK={}
root = 0

def strong_connect_at(s,n,E):
   index = s
   S = []
   indices = {}
   lowlinks = {}
   def strongconnect(index,v):
#       print("strong connect")
#       print(v)
#       print(S)
#       print(indices)
#       print(lowlinks)
       indices[v]=index
       lowlinks[v]=index
       index+=1
       S.append(v)
       for u in E[v]:
           if u not in indices:
               strongconnect(index,u)
               lowlinks[v] = min(lowlinks[v],lowlinks[u])
           elif u in S:
               lowlinks[v] = min(lowlinks[v],indices[u])
       ret = [] 
       if lowlinks[v]==indices[v]:
           while S[-1] != v:
               ret.append(S.pop())
           ret.append(S.pop())
       print("ret")
       print(ret)
#       print(v)
#       print(S)
#       print(indices)
#       print(lowlinks)
       return ret
   return strongconnect(s,s)

def adjacency(n,E):
    adj = {}
    for i in range(n):
        adj[i]=[]
    for (s,t) in E:
        adj[s].append(t)
    return adj

def induced(vs,E):
    rc = {}
    for i in vs:
        rc[i]=[x for x in E[i] if x in vs]
    return rc

def write_cycle(x):
    if len(x)>5:
        #print("cycle")
        path = [reverse_nodes_mapping[i].name for i in x]
        #print(path)
        if 'ADP' not in path and 'ATP' not in path and 'H' not in path and 'H2O' not in path:
            f.write(str(path))
            f.write("\n")

def find_cycles(n,E): # n is the number of vertices in the graph. E is the set of directed edges (i,j) in the graph.
    E = adjacency(n,E)
    for s in range(n):
        Es = induced(range(s,n),E)
        vs = strong_connect_at(s,n,Es)
        Ak = induced(vs,Es)
        B = {}
        blocked = set()
        stack = []
        def unblock(B,v):
            if v in blocked:
                blocked.remove(v)
            if v in B:
                for w in B[v]:
                    B[v].remove(w)
                    if w in blocked:
                        unblock(B,w)

        def circuit(B,v):
            if len(stack)>990:
                print(len(stack))
            found = False
            stack.append(v)
            blocked.add(v)
            for w in Ak[v]:
                if w == s:
                    write_cycle(stack+[s])
                    found = True
                elif w not in blocked:
                    if circuit(B,w):
                        found = True
            if found:
                unblock(B,v)
            else:
                for w in Ak[v]:
                    if w not in B:
                        B[w]=[v]
                    elif v not in B[w]:
                        B[w].append(v)
            stack.pop()
            return found
        circuit(B,s)
reverse_nodes_mapping = dict()
map(lambda k, v: reverse_nodes_mapping.update({k: v}), nodes_mapping.values(), nodes_mapping.keys())

    
print("find cycles:")
f = open("cycles.txt",'w')
find_cycles(n,E)
f.close()
