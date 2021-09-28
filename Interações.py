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
import pandas as pd
import seaborn as sns 
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


def plip(arquivo):
    tree = ET.parse(arquivo)
    root = tree.getroot()
    if len(root.findall('.//bindingsite')) != 1:
        print(arquivo)
        raise Exception("MAIS DE UM BINDING SITE FOI ENCONTRADO!")
        
    props = root.find('.//num_rotatable_bonds').text    
    
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
    
    
caminho = glob('/home/gessualdo/gessualdo/peptbase_1/*/report.xml')
pdb_ids = []
num_rotatable_bonds = []
hydrophobic_interaction = []
hydrogen_bond = []
water_bridge = []
pi_stack = []
salt_bridge = []
pi_cation_interaction = []
halogen_bond = []
metal_complex = []

for report1 in caminho:
    fp1, props = plip(report1)
    num_rotatable_bonds.append(props)
    pdb_ids.append(report1.split('/')[-2])

    hydrophobic_interaction_count = 0
    hydrogen_bond_count = 0
    water_bridge_count = 0
    pi_stack_count = 0
    salt_bridge_count = 0
    pi_cation_interaction_count = 0
    halogen_bond_count = 0
    metal_complex_count = 0
    
    for inter in fp1:
        resnr, reschain, typ = inter
        if typ == 'hydrophobic_interaction':
            hydrophobic_interaction_count = hydrophobic_interaction_count + 1
        elif typ == 'hydrogen_bond':
            hydrogen_bond_count = hydrogen_bond_count + 1
        elif typ == 'water_bridge':
            water_bridge_count = water_bridge_count + 1
        elif typ == 'pi_stack':
            pi_stack_count = pi_stack_count + 1
        elif typ == 'salt_bridge':
            salt_bridge_count = salt_bridge_count + 1
        elif typ == 'pi_cation_interaction':
            pi_cation_interaction_count = pi_cation_interaction_count + 1
        elif typ == 'halogen_bond':
            halogen_bond_count = halogen_bond_count + 1
        elif typ == 'metal_complex':
            metal_complex_count = metal_complex_count + 1
        else:
            raise Exception('TIPO DE INTERAÇÃO DESCONHECIDA')


    hydrophobic_interaction.append(hydrophobic_interaction_count)
    hydrogen_bond.append(hydrogen_bond_count)
    water_bridge.append(water_bridge_count)
    pi_stack.append(pi_stack_count)
    salt_bridge.append(salt_bridge_count)
    pi_cation_interaction.append(pi_cation_interaction_count)
    halogen_bond.append(halogen_bond_count)
    metal_complex.append(metal_complex_count)

data = {
    'PDB.ID':pdb_ids, 
    'num_rotatable_bonds':num_rotatable_bonds,
    'hydrophobic_interaction':hydrophobic_interaction,
    'hydrogen_bond':hydrogen_bond,
    'water_bridge':water_bridge,
    'pi_stack':pi_stack,
    'salt_bridge':salt_bridge,
    'pi_cation_interaction':pi_cation_interaction,
    'halogen_bond':halogen_bond,
    'metal_complex':metal_complex
}

df1 = pd.DataFrame(data=data)


# In[3]:


df1.columns


# In[8]:


df2 = pd.read_csv('/media/gessualdo/C4A81892A8188558/PeptBase.csv')
df2.columns


# In[9]:


peptbase = pd.merge(df1,df2, on = 'PDB.ID')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:



    


# In[ ]:





# In[ ]:




