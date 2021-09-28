#!/usr/bin/env python
# coding: utf-8

# In[4]:


import xml.etree.ElementTree as ET
from xml.dom import minidom
import numpy as np
import seaborn as sb
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd


# In[5]:


def plip(arquivo):
    tree = ET.parse(arquivo)
    root = tree.getroot()
    if len(root.findall('.//bindingsite')) != 1:
        print(arquivo)
        raise Exception("MAIS DE UM BINDING SITE FOI ENCONTRADO!")
        
    fingerprint = set()
    for inter_type in ['hydrophobic_interaction','hydrogen_bond','water_bridge','pi_stack','salt_bridge',
                       'pi_cation_interaction','halogen_bond','metal_complex']:
        for inter in root.findall('.//' + inter_type):
            resnr = inter.find('./resnr').text
            reschain = inter.find('./reschain').text
            typ = inter.tag
            tup = (resnr, reschain, typ)
            fingerprint.add(tup)
    return fingerprint


# In[6]:


caminho = glob('/home/gessualdo/gessualdo/peptbase_1/*/report.xml')


# In[7]:


def tversky(a,b,alfa,beta):
    tversky = len((a).intersection(b))/(len((a).intersection(b))+alfa*len(a-b)+beta*len(b-a))
    
    return (tversky)

rotulos = []
matriz = []
cont_1 = 0
for report1 in caminho:
    fp1 = plip(report1)
    cont_1 = cont_1 +1
    rotulos.append(report1.split('/')[5])
    matriz.append([])
    cont_2 = 0
    for report2 in caminho:    
        fp2 = plip(report2)
        cont_2 = cont_2 +1
        tan = tversky(fp1,fp2, 0.5, 0.5)
        matriz[-1].append(tan)
        if tan > 0.7 and report1 != report2:
            print(report1.split('/')[5],tan,report2.split('/')[5])
print(matriz)
          
         
            
            
        


# In[6]:


#plt.figure(figsize=(23,18))
#plt.title('P24941_series1.Similaridade das interações entre fragmentos e ligantes.' )
#mask = matriz
#sb.heatmap(matriz, xticklabels = rotulos, yticklabels = rotulos, cmap="vlag", linewidths=.5, fmt = 'd')
#plt.show()


# In[ ]:





# In[7]:


#g = sb.clustermap(matriz,
                   #figsize=(18, 15),
                   #row_cluster=False, yticklabels = rotulos, xticklabels = rotulos, cmap = 'vlag',linewidths=.5, dendrogram_ratio=(.1, .2),
                   #cbar_pos=(0, .2, .03, .4))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




