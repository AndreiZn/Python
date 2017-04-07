# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 20:14:28 2017

@author: Andrey
"""

#lesson 3
import pandas as pd

#df = pd.read_csv('ZILL-Z77006_LPC.csv')
#df.set_index ('Date', inplace = True)
#
#df.to_csv('newcsv2.csv')
#
##print (df.head())
#
#df = pd.read_csv('newcsv2.csv', index_col = 'Date')
##print (df.head())
##csv doesn't have index
#
##rename columns:
#df.columns = ['Austin_HPI']
##print (df.head())
#df.to_csv('newcsv3.csv')
##saving without headers 
#df.to_csv('newcsv4.csv', header = False)
#
#df = pd.read_csv('newcsv4.csv', names=['Date', 'Austin_HPI'], index_col = 0)
#print (df.head())

#df.to_html('example.html')

df = pd.read_csv('newcsv4.csv', names=['Date', 'Austin_HPI'])
df.rename (columns={'Austin_HPI':'77006_HPI'}, inplace = True)
print(df.head())