#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 00:02:25 2017

@author: andreyzn
"""

from __future__ import division #so that 3/2 = 1.5, but not 1 as in Python2.* 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


maxcoord = 64
numofcells = 40
hx, hy, hz = maxcoord/numofcells, maxcoord/numofcells, maxcoord/numofcells
h = hx #now hx = hy = hz = h
R = 1.999*hx
V = 4/3*np.pi*R**3 

#reading the halo DataFrame from the file
halo_df = pd.read_csv('halo_df_64Mpc.csv', index_col=0)

#choosing halos with mass>100*prtclmass = 1e9
prtclmass = 2.08108e07
halo_df.sort_values(['M200c'], inplace=True)
halo_df = halo_df.loc[halo_df['M200c'] > 1e3*prtclmass]

#normalized density:
temp_arr = np.array(halo_df.Density)/V
halo_df.loc[:, 'Density'] = pd.Series(temp_arr, index=halo_df.index)

#aver_dens = np.average(halo_df.Density)
aver_dens = 512**3/maxcoord**3
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
list_of_dframes = np.array_split(halo_sorted,   num_of_bins)   

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
    #plot:

    fig1 = plt.figure(1)
    fig1.set_size_inches(15.5, 8.5)
        
    log = 1
    if log==1:
       scale = 'logscale'
       plt.xscale('log')
       plt.yscale('log')
    else:
       scale = 'normal scale' 
       
    ttl = 'Number of halos with mass greater then M. '+str(num_of_bins)+' contrast bins. Cube size='+str(maxcoord)+'Mpc. R='+str(round(R, 1))+'Mpc.'+' Scale='+scale
    plt.title(ttl)
    
    plt.xlabel('M, Msun/h')
    plt.ylabel('n(>M)')    
    plt.plot (M, n)
    leg[i] = 'contr=' + str(round(np.average(df_sorted.Contrast),3))

#legend
plt.legend(leg)
plt.show()






# differentiation 

fig2 = plt.figure(2)
fig2.set_size_inches(15.5, 8.5)

deriv = np.diff(n)/np.diff(np.log(M))
deriv = np.squeeze(deriv) #otherwise its shape is (1,11776)
deriv = np.append(deriv, deriv[-1]) #to make it the same size as M

log = 0
if log==1:
   scale = 'logscale'
   plt.xscale('log')
   plt.yscale('log')
else:
   scale = 'normal scale' 

ttl = 'Mass function. Cube size='+str(maxcoord)+'Mpc. R='+str(round(R, 1))+'Mpc.'+' Scale='+scale
plt.title(ttl)
   
plt.xlabel('M, Msun/h')
plt.ylabel('dn(>M)/dlogM')    

plt.plot (M, deriv)

plt.show()