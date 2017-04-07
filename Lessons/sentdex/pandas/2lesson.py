# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 19:24:25 2017

@author: Andrey
"""

import datetime
import pandas.io.data as web
from matplotlib import style
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

web_stats = {'Day':[1,2,3,4,5,6],
             'Visitors':[43,34,65,56,29,76],
             'Bounce Rate':[65,67,78,65,45,52]}

df = pd.DataFrame(web_stats)

#print(df.head())
#print(df.tail())
#print (df.tail(2))

#df = df.set_index('Day')
#print(df.head())

df.set_index('Day', inplace = True)
#print(df)
#print (df['Visitors'])
#print (df.Visitors)

#print(df[['Bounce Rate', 'Visitors']])

print (np.array(df[['Bounce Rate', 'Visitors']]))

df2 = pd.DataFrame(np.array(df[['Bounce Rate', 'Visitors']]))

print (df2)