#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')

import xml.etree.ElementTree as ET
from xml.dom import minidom
import numpy as np
import seaborn as sb
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage


# In[2]:


def plip(arquivo):
    tree = ET.parse(arquivo)
    root = tree.getroot()
    if len(root.findall('.//bindingsite')) != 1:
        print(arquivo)
        raise Exception("MAIS DE UM BINDING SITE FOI ENCONTRADO!")
     
    num_heavy_atoms = root.find('.//num_heavy_atoms').text
    num_hbd = root.find('.//num_hbd').text
    num_unpaired_hbd = root.find('.//num_unpaired_hbd').text
    num_hba = root.find('.//num_hba').text
    num_unpaired_hba = root.find('.//num_unpaired_hba').text
    num_aromatic_rings = root.find('.//num_aromatic_rings').text
    num_rotatable_bonds = root.find('.//num_rotatable_bonds').text
    molweight = root.find('.//molweight').text
    logp = root.find('.//logp').text
    
    props = {
        'num_heavy_atoms':int(num_heavy_atoms),
        'num_hbd': int(num_hbd),
        'num_unpaired_hbd':int(num_unpaired_hbd),
        'num_hba':int(num_hba),
        'num_unpaired_hba':int(num_unpaired_hba),
        'num_aromatic_rings':int(num_aromatic_rings),
        'num_rotatable_bonds':int(num_rotatable_bonds),
        'molweight':float(molweight),
        'logp':float(logp)    
    }
    
    fingerprint = []
    for inter_type in ['hydrophobic_interaction','hydrogen_bond','water_bridge','pi_stack','salt_bridge',
                       'pi_cation_interaction','halogen_bond','metal_complex']:
        for inter in root.findall('.//' + inter_type):
            resnr = inter.find('./resnr').text
            reschain = inter.find('./reschain').text
            typ = inter.tag
            tup = (resnr, reschain, typ)
            fingerprint.append(tup)
    return fingerprint, props



def tversky(x,y,alfa,beta):
    if x == set() and y == set():
        tversky = 0
    else:
        tversky = len((x).intersection(y))/(len((x).intersection(y))+alfa*len(x-y)+beta*len(y-x))
    
    return (tversky)


# In[50]:


peptbase = pd.read_csv('/media/gessualdo/C4A81892A8188558/PeptBase.csv')

peptbase_group = peptbase.groupby(['UniProt.Accession'])
Uniprot = filtered = peptbase_group.filter(lambda x: x['UniProt.Accession'].count() > 5)


# In[23]:


Uniprot


# In[51]:


protein = str(input('Código Uniprot?'))
pdb = Uniprot[Uniprot['UniProt.Accession'] == protein]
pdbs = list(pdb['PDB.ID'])


# In[57]:



   


# In[88]:


peptbase = pd.read_csv('/media/gessualdo/C4A81892A8188558/PeptBase.csv')

peptbase_group = peptbase.groupby(['UniProt.Accession'])
Uniprot = filtered = peptbase_group.filter(lambda x: x['UniProt.Accession'].count() > 5)




for uniprot in Uniprot['UniProt.Accession'].unique():
    print(uniprot)
    pdb = Uniprot[Uniprot['UniProt.Accession'] == uniprot]
    pdbs = list(pdb['PDB.ID']) 

    count_resi = {}

    rotulos = pdbs
    mat = np.zeros((len(pdbs), len(pdbs)))
    count1 = 0
    for pdb1 in pdbs:
        caminho1 = '/home/gessualdo/gessualdo/peptbase_1/' + pdb1 + '/report.xml'
        count2 = 0
        for pdb2 in pdbs:
            caminho2 = '/home/gessualdo/gessualdo/peptbase_1/' + pdb2 + '/report.xml'
            fp1, _ = plip(caminho1)
            fp2, _ = plip(caminho2)
            fp1 = set(fp1)
            fp2 = set(fp2)
            mat[count1, count2] = 1 - (tversky(min(fp1,fp2),max(fp2,fp1), 0.9, 0.1))
            residuos = set()

            count2 = count2 + 1

        count1 = count1 + 1

        for resid, chain, _ in fp1:
            res = f'{resid}{chain}'
            if res not in count_resi:
                count_resi[res] = +1
            else:
                count_resi[res] += 1

    print(count_resi)
    aminoacidos = list(count_resi.values()) 
    plot = pd.Series(aminoacidos).hist()
    fig = plot.get_figure()
    plt.title(f'Uniprot:{uniprot}')
    plt.xlabel('Frequency of amino acids present in each Uniprot')
    plt.ylabel('Count')
    fig.savefig(f"/home/gessualdo/gessualdo/histograma/{uniprot}.png")
    plt.clf()
    
    
    Z = linkage(mat, 'ward')
    fig = plt.figure(figsize=(10, 10))
    dn = dendrogram(Z,  labels=rotulos)
   
    plt.title(uniprot)
    fig.savefig(f"/home/gessualdo/gessualdo/clustmap/{uniprot}.png")
    plt.clf()
    


        
        
    
  


# In[53]:


pd.Series(aminoacidos).hist();


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[78]:




sb.clustermap(mat,
                   figsize=(10, 10),
                   row_cluster=True, yticklabels = rotulos, xticklabels = rotulos, linewidths=.0, cmap ='Blues', dendrogram_ratio=(.1, .2),
                   cbar_pos=(1, .2, .03, .4), method = 'ward', standard_scale = 1)

                    


# In[ ]:





# In[8]:


#complete, ward, single
Z = linkage(mat, 'ward')
fig = plt.figure(figsize=(10, 10))
dn = dendrogram(Z,  labels=rotulos)
plt.show()


# In[10]:


Z = linkage(mat, 'complete')
fig = plt.figure(figsize=(10, 10))
dn = dendrogram(Z,  labels=rotulos)
plt.show()


# In[9]:


Z = linkage(mat, 'single')
fig = plt.figure(figsize=(10, 10))
dn = dendrogram(Z,  labels=rotulos)
plt.show()


# In[ ]:




