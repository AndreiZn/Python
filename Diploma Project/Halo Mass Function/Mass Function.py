#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 21:28:53 2017

@author: andreyzn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#reading the halo DataFrame from the file
halo_df = pd.read_csv('./32Mpc/halo_df.csv', index_col=0)

#choosing halos with mass>100*prtclmass = 1e9
prtclmass = 2.08108e07
halo_df.sort_values(['M200c'], inplace=True)
halo_df = halo_df.loc[halo_df['M200c'] > 1e3*prtclmass]

#
aver_dens = np.average(halo_df.Density)
halo_df.loc[:, 'Contrast'] = (halo_df.Density - aver_dens)/aver_dens

#sorting the halo DataFrame
halo_sorted = halo_df.sort_values(['Contrast'])

#Index is mixed up. Changing it to 0, 1, 2, ... int_array
int_array = np.array([i for i in range(np.size(halo_sorted, 0)) ])
halo_sorted.loc[:,'new_index'] = pd.Series(int_array, index = halo_sorted.index)
halo_sorted.set_index('new_index', inplace=True)
#halo_sorted.set_index(pd.Series(int_array, index = halo_sorted.index), inplace=True)


aver_contrast = np.average(halo_sorted.Contrast)
sigma2 = np.average((halo_sorted.Contrast- aver_contrast)**2/aver_contrast)
#sigma = np.sqrt(sigma2)

halosize = np.size(halo_sorted, 0)

num_of_bins = 4

#dividing halos into bins
list_of_dframes = np.array_split(halo_sorted, num_of_bins)   

#fig = plt.figure()

#leg is used in legend
leg = ['']*num_of_bins

for i in range (num_of_bins):
    df = list_of_dframes[i]
    df_sorted = df.sort_values(['M200c'])
    size = np.size(df_sorted,0) 
    int_array = np.array([j for j in range(0, size) ])   
    int_array = int_array[::-1]
    df_sorted.set_index(pd.Series(int_array, index = df_sorted.index), inplace=True)
    
    #values to plot:
    n, M = np.array(df_sorted.index), np.array(df_sorted.M200c)
    plt.title('Mass function for different bins')
    plt.xlabel('log M, Msun/h')
    plt.ylabel('log n(>M)')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot (M, n)
    leg[i] = 'contr=' + str(round(np.average(df_sorted.Contrast),3))

#legend
plt.legend(leg)
plt.show()