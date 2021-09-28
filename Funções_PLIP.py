#!/usr/bin/env python
# coding: utf-8

# **PeptBase: Proteínas-peptídeos**

# In[1]:


#Importando bibliotecas

import xml.etree.ElementTree as ET
from xml.dom import minidom
import numpy as np
import seaborn as sb
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import Counter
        


# **Funções para obtenção das propriedades Físico-químicas da Peptbase**

# In[ ]:





# In[76]:


#obter as características físico-químicas dos peptídeos.


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



# In[ ]:


#Métrica que avalia a similaridade dos resíduos de aminoácidos.


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
        #if tan > 0.7 and report1 != report2:
            #print(report1.split('/')[5],tan,report2.split('/')[5])


          


# In[10]:


#Verificar o tamanho de todos os peptídeos, espera-se que todos sejam abaixo de 10 aa.  

def size_peptide(arquivo):
    tree = ET.parse(arquivo)
    root = tree.getroot()
    if len(root.findall('.//bindingsite')) != 1:
        raise Exception("MAIS DE UM BINDING SITE FOI ENCONTRADO!")
  
    longname = root.find('.//longname').text
    
    return len(str(longname).replace("-",',').split(','))


# In[3]:


#Pegar todos os resíduos de aminoácidos dos peptídeos

def aminoacidos(arquivo): 
    tree = ET.parse(arquivo)
    root = tree.getroot()
    if len(root.findall('.//bindingsite')) != 1:
        raise Exception("MAIS DE UM BINDING SITE FOI ENCONTRADO!")
        
    aa = root.find('.//longname').text
        
    return aa.replace('-', ',').split(',')
    


# In[4]:


#Converter todas as unidades de afinidades para nM.


def converter(caminho):
    value_seq = {}
    df = pd.read_csv(caminho)
    for _, row in df.iterrows():
        value = row['Affinity.Value']
        unit = row['Affinity.Original.Unit']
        pdbid = row['PDB.ID']
        
        if unit == 'uM':
            value = value*1000
            
        elif unit == 'nM':
            value = value*1
            
        elif unit == 'mM':
            value = value* 1e+6
            
        elif unit == 'M':
            value = value* 1e+9
            
        elif unit == 'M^-1':
            value = value* 0.1
            
        elif unit == 'pM':
            value = value/1000
        
        else:
            raise Exception ("UNIDADE DESCONHECIDA" + unit)
            
        seq = row['Peptide.Sequence..3.letter.']
        
        value_seq[pdbid] = value, seq.split(' ')
        
    return value_seq
        


# **Análises feitas na PeptBase**

# Verificar se existe algum PDB com peptídeos com mais de 10 aa.

# In[122]:



database = pd.read_csv('/home/gessualdo/Documentos/peptbase_/PeptBase.csv')
pdb_id = df['PDB.ID']
id_pdb = list(dict.fromkeys(pdb_id))
peptide_resid = database['Peptide.Resid']
chain = database['Peptide.Chain']

ID_PDB = []
ID_RES = []
CHAIN_PEPTIDEO = []

for pdb, res, chain_pep in zip(data, peptide_resid, chain):
    caminho1 = '/home/gessualdo/Documentos/resultados plip/' + pdb + '/report.xml'
    if size_peptide(caminho1) > 10:
        ID_PDB.append(pdb)
        ID_RES.append(res)
        CHAIN_PEPTIDEO.append(chain_pep)
        print('PDB.ID:',pdb, 'Tamanho do Peptídeo:', size_peptide(caminho1),'ID.RESID:' ,res)


# Os peptídeos foram oriundos de uma base de dados de complexo de proteínas-petídeos com o máximo de 10 aa. Dessa forma, esses valores estão incorretos, sendo necessários realizar alguns ajustes no PyMOL.

# In[64]:


from pymol import cmd as pm
import subprocess


# In[63]:


CHAIN_PEPTIDEO


# In[65]:


for id_names, resi, chain in zip(ID_PDB, ID_RES, CHAIN_PEPTIDEO):


    resi_number = {
        id_names:(resi ,chain),

    }

    for obj, (resi, chain) in resi_number.items():
        pm.fetch(obj)
        pm.alter(f'bymolecule {obj} and resi {resi} and chain {chain}', 'chain="P"')
        arquivo_pdb = f'/home/gessualdo/Documentos/resultados plip/{obj}/{obj}.pdb'
        pm.save(arquivo_pdb, obj)

        print(obj, chain)
        caminho_plip = '/home/gessualdo/anaconda3/envs/envname/bin/plip'
        plip = subprocess.run([caminho_plip, '-f', arquivo_pdb, '--peptide', 'P',
                               '-o', f'/home/gessualdo/Documentos/resultados plip/{obj}', '-x'])


# In[ ]:





# Obter as propriedas físico-químicas dos peptídeos

# In[123]:


#Endereço com os arquivos XML
caminho = glob('/home/gessualdo/Documentos/resultados plip/*/report.xml')

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

peptbase = {
    'PDB.ID':pdb_ids, 
    'Num.Rotatable.Bonds':num_rotatable_bonds,
    'Hydrophobic.Interaction':hydrophobic_interaction,
    'Hydrogen.Bond':hydrogen_bond,
    'Water.Bridge':water_bridge,
    'Pi.Stack':pi_stack,
    'Salt.Bridge':salt_bridge,
    'Pi.Cation.Interaction':pi_cation_interaction,
    'Halogen.Bond':halogen_bond,
    'Metal.Complex':metal_complex
}

df = pd.DataFrame(data=peptbase)


# In[124]:


database = database.drop(['Unnamed: 0.1'], axis = 1)
database = database.drop(['Unnamed: 0'], axis = 1)


# In[125]:


#Merge as informações com as propriedades físico químicas com o dataset original

peptbase = database.merge(df, left_on='PDB.ID', right_on='PDB.ID')


# In[ ]:





# In[127]:


peptbase.to_csv('/home/gessualdo/Documentos/peptbase_/PeptBase_.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




