#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[9]:


file = pd.read_csv('NEW_SPINS_output_file.csv')
plsc_brain = pd.read_csv('NEW_rs_conn_LV1_correct_with_zeros.csv')

#extracting the first column
column = file.iloc[:, 1]

#multiplying values by 10000
big_column = column*10000

#changes the values all to 0
plsc_brain['r'] = 0 
#attach it to the end of the other csv
plsc_brain['r'] += big_column

#delete the first two columns:
plsc_brain_drop = plsc_brain.drop(plsc_brain.columns[0], axis=1)
#plsc_brain_drop2 = plsc_brain_drop.drop(plsc_brain_drop.columns[0], axis=1) #weird cause it doesnt consider the first unnamed column in the csv to be anything


# In[11]:


pd.DataFrame(plsc_brain_drop).to_csv("NEW_rs_conn_LV1.csv")

