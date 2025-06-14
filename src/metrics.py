# Modified from: original_file.py from Project Name (https://github.com/user/repo)
# Original License: Apache License 2.0
# Changes made by Susanna Grigson on 8 October 2024: Copy functions for metrics: read_pdb, calc_pdockq, calc_pdockq_multichain 
#!/usr/bin/env python 

"""
Methods for computing metrics associated with protein interactions 
"""

# imports for reading in 
import pandas as pd 
import json
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from collections import defaultdict


def parse_atm_record(line):
    '''
    Extract details from line in pdb file 
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['B'] = float(line[60:66])

    return record


def read_pdb(pdbfile):
    '''Read a pdb file predicted with AF and rewritten to conatin all chains
    
    from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
    '''

    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            
            #Get CB - CA for GLY
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                    chain_plddt[record['chain']] = [record['B']]


    #Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt


def calc_pdockq(chain_coords, chain_plddt, t):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    
    from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py?ref_type=heads
    
    param: t - distance in angstroms - defulat to use is 8 
    '''

    #Get coords and plddt per chain
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    #Calc 2-norm
    mat = np.append(coords1, coords2,axis=0)
    a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1,l1:] #upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists<=t)

    if contacts.shape[0]<1:
        pdockq=0
        ppv=0
    else:
        #Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
        #Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt*np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

    return pdockq

def calc_pdockq_multichain(chain_coords, chain_plddt, t = 8):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    
    from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py?ref_type=heads
    
    param: t - distance in angstroms - defulat to use is 8 
    '''

    #Get coords and plddt per chain
    ch_list = [*chain_coords.keys()]
    coords_list = [chain_coords[ch] for ch in ch_list] 
    plddt_list = [chain_plddt[p] for p in ch_list] 
    
    
    # perform pairwise calculations to get all of the pdockq scores 
    pdockq_scores = [] 
    
    for i in range(len(ch_list)): 
        
        for j in range(i + 1,len(ch_list)):

            #Calc 2-norm
            mat = np.append(coords_list[i], coords_list[j],axis=0)
            a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
            dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
            l1 = len(coords_list[i])
            contact_dists = dists[:l1,l1:] #upper triangular --> first dim = chain 1
            contacts = np.argwhere(contact_dists<=t)

            if contacts.shape[0]<1:
                pdockq=0
                ppv=0
            else:
                #Get the average interface plDDT
                avg_if_plddt = np.average(np.concatenate([plddt_list[i][np.unique(contacts[:,0])], plddt_list[j][np.unique(contacts[:,1])]]))
                #Get the number of interface contacts
                n_if_contacts = contacts.shape[0]
                x = avg_if_plddt*np.log10(n_if_contacts)
                pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018
                
                
            pdockq_scores.append(pdockq)

    return pdockq_scores

#### pdockq2 score 
### This is the code from https://gitlab.com/ElofssonLab/afm-benchmark/-/blob/main/src/pdockq2.py?ref_type=heads 
### have edited it to read a json file rather than a pickle score as input. Gives on score for each chain  

def retrieve_IFplddt(structure, chain1, chain2_lst, max_dist):
    ## generate a dict to save IF_res_id
    chain_lst = list(chain1) + chain2_lst

    ifplddt = []
    contact_chain_lst = []
    for res1 in structure[0][chain1]:
        for chain2 in chain2_lst:
            count = 0
            for res2 in structure[0][chain2]:
                if res1.has_id('CA') and res2.has_id('CA'):
                   dis = abs(res1['CA']-res2['CA'])
                   ## add criteria to filter out disorder res
                   if dis <= max_dist:
                      ifplddt.append(res1['CA'].get_bfactor())
                      count += 1

                elif res1.has_id('CB') and res2.has_id('CB'):
                   dis = abs(res1['CB']-res2['CB'])
                   if dis <= max_dist:
                      ifplddt.append(res1['CB'].get_bfactor())
                      count += 1
            if count > 0:
              contact_chain_lst.append(chain2)
    contact_chain_lst = sorted(list(set(contact_chain_lst)))   


    if len(ifplddt)>0:
       IF_plddt_avg = np.mean(ifplddt)
    else:
       IF_plddt_avg = 0

    return IF_plddt_avg, contact_chain_lst


def retrieve_IFPAEinter(structure, paeMat, contact_lst, max_dist):
    ## contact_lst:the chain list that have an interface with each chain. For eg, a tetramer with A,B,C,D chains and A/B A/C B/D C/D interfaces,
    ##             contact_lst would be [['B','C'],['A','D'],['A','D'],['B','C']]

 
    chain_lst = [x.id for x in structure[0]]
    seqlen = [len(x) for x in structure[0]]
    ifch1_col=[]
    ifch2_col=[]
    ch1_lst=[]
    ch2_lst=[]
    ifpae_avg = []
    d=10
    for ch1_idx in range(len(chain_lst)):
      ## extract x axis range from the PAE matrix
      idx = chain_lst.index(chain_lst[ch1_idx])
      ch1_sta=sum(seqlen[:idx])
      ch1_end=ch1_sta+seqlen[idx]
      ifpae_col = []   
      ## for each chain that shares an interface with chain1, retrieve the PAE matrix for the specific part.
      for contact_ch in contact_lst[ch1_idx]:
        index = chain_lst.index(contact_ch)
        ch_sta = sum(seqlen[:index])
        ch_end = ch_sta+seqlen[index]
        remain_paeMatrix = paeMat[ch1_sta:ch1_end,ch_sta:ch_end]
        #print(contact_ch, ch1_sta, ch1_end, ch_sta, ch_end)        

        ## get avg PAE values for the interfaces for chain 1
        mat_x = -1
        for res1 in structure[0][chain_lst[ch1_idx]]:
          mat_x += 1
          mat_y = -1
          for res2 in structure[0][contact_ch]:
              mat_y+=1
              if res1['CA'] - res2['CA'] <=max_dist:
                 ifpae_col.append(remain_paeMatrix[mat_x,mat_y])
                    
      ## normalize by d(10A) first and then get the average
      if not ifpae_col:
        ifpae_avg.append(0)
      else:
        norm_if_interpae=np.mean(1/(1+(np.array(ifpae_col)/d)**2))
        ifpae_avg.append(norm_if_interpae)

    return ifpae_avg
    

def calc_pmidockq(ifpae_norm, ifplddt):
    df = pd.DataFrame()
    df['ifpae_norm'] = ifpae_norm
    df['ifplddt'] = ifplddt
    df['prot'] = df.ifpae_norm*df.ifplddt
    fitpopt = [1.31034849e+00, 8.47326239e+01, 7.47157696e-02, 5.01886443e-03] ## from orignal fit function  
    df['pmidockq'] = sigmoid(df.prot.values, *fitpopt)

    return df

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return (y)

def fit_newscore(df, column):

    testdf = df[df[column]>0]

    colval = testdf[column].values
    dockq = testdf.DockQ.values
    xdata =colval[np.argsort(colval)]
    ydata = dockq[np.argsort(dockq)]

    p0 = [max(ydata), np.median(xdata),1,min(ydata)] # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata,p0)# method='dogbox', maxfev=50000)
    
    return popt


def calc_pdockq2(json_scores,pdb, dist = 8):
    """ 
    dist refers to distance used in computation of dockq score. Defualt = 8
    """ 

    pdbp = PDBParser(QUIET=True)

    structure = pdbp.get_structure('', pdb)
    chains = []
    for chain in structure[0]:
        chains.append(chain.id)

    remain_contact_lst=[]
    
    ## retrieve interface plDDT at chain-level
    plddt_lst = []
    for idx in range(len(chains)):
        chain2_lst = list(set(chains)-set(chains[idx]))
        IF_plddt, contact_lst = retrieve_IFplddt(structure, chains[idx], chain2_lst, dist)
        plddt_lst.append(IF_plddt)
        remain_contact_lst.append(contact_lst)

    # Read the JSON file and load its contents into a dictionary
    with open(json_scores, 'r') as json_file:
        data = json.load(json_file)

    avgif_pae = retrieve_IFPAEinter(structure, np.array(data['pae']), remain_contact_lst, dist)
    ## calculate pmiDockQ

    res = calc_pmidockq(avgif_pae, plddt_lst)
    
    return res['pmidockq'], avgif_pae, plddt_lst

def calc_interface_contacts(chain_coords, t = 8):
    '''Calculate the number of interface contacts
    
    param: t - distance in angstroms - default to use is 8 
    '''

    # Get coords and plddt per chain
    ch_list = [*chain_coords.keys()]
    coords_list = [chain_coords[ch] for ch in ch_list] 
    
    # List to store the number of interface contacts
    interface_contacts = [] 
    
    for i in range(len(ch_list)): 
        for j in range(i + 1, len(ch_list)):
            # Calculate 2-norm distances
            mat = np.append(coords_list[i], coords_list[j], axis=0)
            a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
            dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
            l1 = len(coords_list[i])
            contact_dists = dists[:l1, l1:] # Upper triangular --> first dim = chain 1
            contacts = np.argwhere(contact_dists <= t)

            # Get the number of interface contacts
            n_if_contacts = contacts.shape[0]
            interface_contacts.append(n_if_contacts)

    return interface_contacts
