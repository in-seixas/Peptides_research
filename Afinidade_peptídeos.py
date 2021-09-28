#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params , palette = 'rocket')


# In[3]:


df = pd.read_csv('/home/gessualdo/Documentos/peptbase_/PeptBase_.csv')


# In[4]:


df


# In[5]:


df['Pki'] = -np.log10(df['Affinity.Value'])


# In[ ]:





# In[6]:


df.loc[df['Pki'] < df['Pki'].quantile([0.25])[0.25], ['Affinity.Class1']] = 'low'
df.loc[df['Pki'] > df['Pki'].quantile([0.75])[0.75], ['Affinity.Class1']] = 'high'


# In[7]:


df


# In[7]:


sns.boxplot(x = 'Pki', y= 'Affinity.Class1', data = df)


# In[8]:




l = [{},{},{},{},{},{},{},{},{}]
p = [0,0,0,0,0,0,0,0,0]

for peptideo in df['Peptide.Sequence']:
    pep_len = len(peptideo)
    d = l[pep_len - 2]
    for aa in peptideo:
        if aa in d:
            
            d[aa] = d[aa] + 1
            
        else:
            d[aa] = 1
            
for i,d in enumerate(l):
    p[i] = sum(d.values())
    
for i, d in enumerate(l):
    for aa in d:
        d[aa] = (d[aa]/p[i])*100
        


# In[9]:


df_freq = pd.DataFrame(l)

df_freq[pd.isnull(df_freq)] = 0

df_freq['Peptide.Length'] = df_freq.index + 2
df_freq


df_freq = df_freq.melt(id_vars=["Peptide.Length"], 
        var_name="Aminoacid", 
        value_name="Percent")

polaridade = {
    "polar_charged+":['R','H', 'K'],
    "polar_charged-":['D','E'],
    "polar_uncharged":['S', 'T', 'N', 'Q'],
    "special_cases":['C', 'U', 'G', 'P'],
    "hydrophobic":['A', 'I', 'L', 'M', 'F','W', 'Y', 'V'] 
}

for polarity, aminoacids in polaridade.items():
    for aa in aminoacids:
        
        df_freq.loc[df_freq['Aminoacid'] == aa, "Polarity"] = polarity


# In[13]:



df_freq['Polarity'].value_counts()


# In[14]:


g = sns.FacetGrid(df_freq, col="Peptide.Length", height=5, aspect=.7, col_wrap=3)
g.map(sns.barplot, "Aminoacid", "Percent", order='Polarity')


# In[15]:



plt.figure(figsize=(10,6))
a = (pd.Series(l[8]).sort_values()*100)/sum(d.values())


a.plot(kind='barh')


# In[20]:


sns.scatterplot(x= 'Pki', y='Peptide.Length', data=df)


# In[51]:


df_freq.to_csv('/home/gessualdo/Documentos/datasets/df_freq.csv')


# In[ ]:





# In[52]:


df.to_csv('/home/gessualdo/Documentos/datasets/PeptBase.csv')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




