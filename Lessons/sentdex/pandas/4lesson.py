# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 21:14:13 2017

@author: Andrey
"""

#lesson 4
import quandl
import pandas as pd

api_key = 'yjxUqzGvctxnxz16ErHx'

df = quandl.get ('FMAC/HPI_AK', authtoken = api_key)

#print(df.head())
#it's going to make list of dataframes
fifty_states = pd.read_html('https://simple.wikipedia.org/wiki/List_of_U.S._states')

#this is a list:
#print (fifty_states)

#this is a dataframe:
#print (fifty_states[0])

#this is a column:
print (fifty_states[0][0])

for abbv in fifty_states[0][0][1:]:
    print ("FMAC/HPI_"+str(abbv))
    