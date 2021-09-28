#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xml.etree.ElementTree as ET
from xml.dom import minidom
import numpy as np
import seaborn as sb
import os
from glob import glob
import matplotlib.pyplot as plt


# In[ ]:



    


# In[13]:


#Função que abre os arquivos XML e realiza a busca pelas interações para cada ligante
def plip(diretorio):
    tree = ET.parse(diretorio)
    root = tree.getroot()         
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
    


# In[14]:


caminho = glob('/home/gessualdo/gessualdo/peptbase_1/*/report.xml')


# In[15]:


def tanimoto(a, b):
    tanimoto = len((a).intersection(b))/len((a).union(b))
    return (tanimoto)

rotulos = []
matriz = []
cont_1 = 0
for report1 in caminho:
    fp1 = plip(report1)
    cont_1 = cont_1 +1
    rotulos.append(report1.split('/')[6])
    matriz.append([])
    cont_2 = 0
    for report2 in caminho:    
        fp2 = plip(report2)
        cont_2 = cont_2 +1
        tan = tanimoto(fp1,fp2)
        matriz[-1].append(tan)
        if tan >= 0.6 and report1 != report2:
            print(report1.split('/')[6],tan,report2.split('/')[6])
            
        
        
        


    
    


# In[ ]:


plt.figure(figsize=(23,18))
plt.title('Q13526_series2.Similaridade das interações entre fragmentos e ligantes.' )
mask = matriz
sb.heatmap(matriz, xticklabels = rotulos, yticklabels = rotulos, cmap="vlag", linewidths=.5, fmt = 'd')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[71]:





# In[19]:



    


# In[22]:





# 

# In[ ]:



    


# In[ ]:





# In[ ]:





# In[37]:





# In[47]:



    



# In[ ]:



        


# In[ ]:





# In[ ]:




