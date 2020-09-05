import numpy as np
import matplotlib.pyplot as plt
from FEM_beams import *


count_node = 9
count_elem = count_node-1
length = 1.0

list_b=[0.0,0.0191,0.0354,0.0461,0.05,0.0461,0.0354,0.0191]
list_k=[0.0191,0.0163,0.0107,0.0039,-0.0039,-0.0107,-0.0163,-0.0191]

list_b_big=np.zeros(30)
list_k_big=np.zeros(30)

for i in range(0,20):
    list_b_big[i]=0.05
    list_k_big[i]=0

dict_nodes={} #nid:[x,y]
dict_elems={} #eid:[type,nid1,nid2,pid]
dict_params={} #pid:{all the required parameters for the element}
for i in range(0,count_node):
    dict_nodes[i]=[length/count_elem*i,0]
for i in range(0,count_elem):
    dict_elems[i]=['BEAM',i,i+1,1]
dict_params[1] = {'cs-profile':'single-layer','E':8.3467e10,'DENS':7.6e3,'B':5e-2,'H':1.5e-1}
list_bc = ['0UX','0UY','%dUX'%(count_node-1,),'%dUY'%(count_node-1,)]

FEM_K,FEM_M,DOFmappingDic =  Assembling_PZTSC(dict_elems,dict_nodes,dict_params,list_b,list_k)

Apply_mecha_BC(FEM_K,FEM_M,DOFmappingDic,list_bc)
mfreqs, mshapes= Solve_MODAL(FEM_K,FEM_M, EigCount=8)
print(mfreqs)



mid = 1
PlotDef(dict_nodes,mshapes[:,mid],DOFmappingDic,title_text='Modal shape (No.%d)'%(mid,))
plt.show()