import matplotlib.pyplot as plt
import numpy as np
fy_paper_label_font = {'family': 'serif','color':  'darkred','weight': 'normal','fontsize': 16}
from matplotlib.ticker import FormatStrFormatter

def GetElementDOF(etype):
    '''
    return a list of DOF labels for each node of the give element type [etype]
    '''
    ElementDOFDic = {'BEAM': ('UX', 'UY', 'THXY'),
                     'PZTBEAM': ('UX', 'UY', 'THXY', 'VT', 'VB')}
    return ElementDOFDic[etype]

def DistributeDOF(dict_elems):
    DOFmappingDic = {}
    eqn = 0
    donelist = []
    for elem in dict_elems.values():
        DOFList = GetElementDOF(elem[0])
        for nID in [elem[1],elem[2]]:
            for DOF in DOFList:
                tag = str(nID) + DOF
                if tag not in donelist:
                    DOFmappingDic[tag] = eqn
                    eqn += 1
                    donelist.append(tag)
    DOFCount = len(donelist)
    return DOFCount,DOFmappingDic


def Assembling(dict_elems,dict_nodes,dict_params):
    DOFCount, DOFmappingDic = DistributeDOF(dict_elems)
    FEM_K = np.zeros((DOFCount, DOFCount), dtype='float')
    FEM_M = np.zeros_like(FEM_K)
    for elem in dict_elems.values():
        # collecting parameter to generate element matrix
        params = dict_params[elem[3]]
        L = np.abs(dict_nodes[elem[1]][0] - dict_nodes[elem[2]][0])
        # obtain element matrix
        if elem[0] == 'BEAM':
            eK = Beam_GetK(params,L)
            eM = Beam_GetM(params,L)
            str0,str1 = str(elem[1]),str(elem[2])
            eDOF = (str0 + 'UX', str1 + 'UX', str0 + 'UY', str0 + 'THXY', str1 + 'UY', str1 + 'THXY')
        elif elem[0] == 'PZTBEAM':
            pass
        # obtain corresponding equation numbers for each element DOF
        eDOFCount = len(eDOF)
        eEqlist = []
        for DOF in eDOF:
            eEqlist.append(DOFmappingDic[DOF])
        # assembling
        for ei in range(0, eDOFCount):
            for ej in range(0, eDOFCount):
                oi = eEqlist[ei]
                oj = eEqlist[ej]
                FEM_K[oi, oj] += eK[ei, ej]
                FEM_M[oi, oj] += eM[ei, ej]
    return FEM_K,FEM_M,DOFmappingDic

def Assembling_PZTOC(dict_elems,dict_nodes,dict_params,list_b,list_k):
    DOFCount, DOFmappingDic = DistributeDOF(dict_elems)
    FEM_K = np.zeros((DOFCount+1, DOFCount+1), dtype='float')
    FEM_M = np.zeros_like(FEM_K)
    
    elemCount_ori=0
    for elem in dict_elems.values():
        elemCount=elemCount_ori*3
        # collecting parameter to generate element matrix
        params = dict_params[elem[3]]
        L = np.abs(dict_nodes[elem[1]][0] - dict_nodes[elem[2]][0])
        # obtain element matrix
        if elem[0] == 'BEAM':
            eK = Beam_GetK(params,L)
            eM = Beam_GetM(params,L)
            str0,str1 = str(elem[1]),str(elem[2])
            eDOF = (str0 + 'UX', str1 + 'UX', str0 + 'UY', str0 + 'THXY', str1 + 'UY', str1 + 'THXY')
        elif elem[0] == 'PZTBEAM':
            pass
        # add pzt parameters to matrix
        U=9.3*0.05/L
        C=9.558e-9
        
        FEM_K[elemCount+1,DOFCount]+=2*U*list_k[elemCount_ori]
        FEM_K[elemCount+2,DOFCount]+=-2*U*L*list_b[elemCount_ori]
        FEM_K[elemCount+4,DOFCount]+=-2*U*list_k[elemCount_ori]
        FEM_K[elemCount+5,DOFCount]+=2*U*L*(list_b[elemCount_ori]+list_k[elemCount_ori])
        
        FEM_K[DOFCount,elemCount+1]+=2*U*list_k[elemCount_ori]
        FEM_K[DOFCount,elemCount+2]+=-2*U*L*list_b[elemCount_ori]
        FEM_K[DOFCount,elemCount+4]+=-2*U*list_k[elemCount_ori]
        FEM_K[DOFCount,elemCount+5]+=2*U*L*(list_b[elemCount_ori]+list_k[elemCount_ori])
        
        FEM_K[DOFCount,DOFCount]+=C*(2*list_b[elemCount_ori]+list_k[elemCount_ori])
        # obtain corresponding equation numbers for each element DOF
        eDOFCount = len(eDOF)
        eEqlist = []
        for DOF in eDOF:
            eEqlist.append(DOFmappingDic[DOF])
        # assembling
        for ei in range(0, eDOFCount):
            for ej in range(0, eDOFCount):
                oi = eEqlist[ei]
                oj = eEqlist[ej]
                FEM_K[oi, oj] += eK[ei, ej]
                FEM_M[oi, oj] += eM[ei, ej]
        elemCount_ori+=1
                
    return FEM_K,FEM_M,DOFmappingDic

def Assembling_PZTSC(dict_elems,dict_nodes,dict_params,list_b,list_k):
    DOFCount, DOFmappingDic = DistributeDOF(dict_elems)
    FEM_K = np.zeros((DOFCount+1, DOFCount+1), dtype='float')
    FEM_M = np.zeros_like(FEM_K)
    
    FEM_K[DOFCount,DOFCount]+=1
    for elem in dict_elems.values():
        # collecting parameter to generate element matrix
        params = dict_params[elem[3]]
        L = np.abs(dict_nodes[elem[1]][0] - dict_nodes[elem[2]][0])
        # obtain element matrix
        if elem[0] == 'BEAM':
            eK = Beam_GetK(params,L)
            eM = Beam_GetM(params,L)
            str0,str1 = str(elem[1]),str(elem[2])
            eDOF = (str0 + 'UX', str1 + 'UX', str0 + 'UY', str0 + 'THXY', str1 + 'UY', str1 + 'THXY')
        elif elem[0] == 'PZTBEAM':
            pass
        # obtain corresponding equation numbers for each element DOF
        eDOFCount = len(eDOF)
        eEqlist = []
        for DOF in eDOF:
            eEqlist.append(DOFmappingDic[DOF])
        # assembling
        for ei in range(0, eDOFCount):
            for ej in range(0, eDOFCount):
                oi = eEqlist[ei]
                oj = eEqlist[ej]
                FEM_K[oi, oj] += eK[ei, ej]
                FEM_M[oi, oj] += eM[ei, ej]
    return FEM_K,FEM_M,DOFmappingDic

def Apply_mecha_BC(FEM_K,FEM_M,DOFmappingDic,list_bc):
    for bc_dof in list_bc:
        eqn = DOFmappingDic[bc_dof]
        temp = FEM_K[eqn, eqn]
        FEM_K[eqn, :] = 0
        FEM_K[:, eqn] = 0
        FEM_M[:, eqn] = 0
        FEM_M[eqn, :] = 0
        FEM_K[eqn, eqn] = temp


###########################################
#  starting - Pure mechanical beams
###########################################
def Beam_GetM(params,l):
    if params['cs-profile'] == 'single-layer':
        dens,b,h = params['DENS'],params['B'],params['H']
        m = dens * b * h * l
    eM = np.array([[1.0 / 3, 1.0 / 6, 0, 0, 0, 0],
                   [1.0 / 6, 1.0 / 3, 0, 0, 0, 0],
                   [0, 0, 156, 22 * l, 54, -13 * l],
                   [0, 0, 22 * l, 4 * l * l, 13 * l, -3 * l * l],
                   [0, 0, 54, 13 * l, 156, -22 * l],
                   [0, 0, -13 * l, -3 * l * l, -22 * l, 4 * l * l]], dtype='float')
    eM[2:, 2:] = eM[2:, 2:] / 420.0
    eM[:, :] = m * eM[:, :]
    return eM

def Beam_GetK(params,l):
    if params['cs-profile'] == 'single-layer':
        E,b,h = params['E'],params['B'],params['H']
        GA = b * h
        GI = b * np.power(h, 3)/12
        kl = E * GA / l
        kb = E * GI/ np.power(l, 3)
    eK = np.array([[1, -1, 0, 0, 0, 0],
                   [-1, 1, 0, 0, 0, 0],
                   [0, 0, 12, 6 * l, -12, 6 * l],
                   [0, 0, 6 * l, 4 * l * l, -6 * l, 2 * l * l],
                   [0, 0, -12, -6 * l, 12, -6 * l],
                   [0, 0, 6 * l, 2 * l * l, -6 * l, 4 * l * l]], dtype='float')
    eK[0:2, 0:2] = kl * eK[0:2, 0:2]
    eK[2:, 2:] = kb * eK[2:, 2:]
    return eK
###########################################
#  end - Pure mechanical beams
###########################################

            
def Solve_MODAL(FEM_K,FEM_M, EigCount):
    import scipy.sparse.linalg as slg
    evals, evecs = slg.eigsh(FEM_K,k=EigCount,M=FEM_M,sigma = 0.0)
    mfreqs = np.power(np.real(evals),0.5)/np.pi/2
    return mfreqs,np.real(evecs)               
            


def PlotDef(dict_nodes,deform,DOFmappingdic,title_text=''):
    xlist = []
    ylist = []
    thxylist = []
    uxlist = []
    for nID, cords in dict_nodes.items():
        xlist.append(cords[0])
        ylist.append(deform[DOFmappingdic[str(nID)+'UY']])
        thxylist.append(deform[DOFmappingdic[str(nID)+'THXY']])
        uxlist.append(deform[DOFmappingdic[str(nID)+'UX']])
    
    arg = np.argsort(xlist)
    xlist = np.array(xlist)[arg]
    ylist = np.array(ylist)[arg]
    thxylist = np.array(thxylist)[arg]
    uxlist = np.array(uxlist)[arg]
    
    plt.figure()
    plt.suptitle(title_text,fontdict = fy_paper_label_font)
    plt.subplot(211)
    plt.ylabel('UY',fontdict = fy_paper_label_font)
    plt.plot(xlist,ylist,lw=1,marker='s',color='green')
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%3.0e'))
    plt.subplot(212)
    plt.plot(xlist,thxylist,lw=1,marker='s',color='orange')
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%3.0e'))
    plt.xlabel('x coordinate /m',fontdict = fy_paper_label_font)
    plt.ylabel('THXY',fontdict = fy_paper_label_font)
    plt.tight_layout()

    plt.figure()
    plt.suptitle(title_text,fontdict = fy_paper_label_font)
    plt.plot(xlist,uxlist,lw=1,marker='s',color='purple')
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%3.0e'))
    plt.xlabel('x coordinate /m',fontdict = fy_paper_label_font)
    plt.ylabel('UX',fontdict = fy_paper_label_font)              
    plt.tight_layout()
                
            
                
            
            
    
            
            
    