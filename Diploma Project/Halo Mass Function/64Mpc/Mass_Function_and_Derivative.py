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
halo_df = pd.read_csv('./Halo dataframes for different radii/halo_df_64Mpc_R_3.2.csv', index_col=0)

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

num_of_bins = 1

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
    N, M = np.array(df_sorted.index), np.array(df_sorted.M200c)
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
    plt.ylabel('N(>M)')    
    plt.plot (M, N, '.')
    leg[i] = 'contr=' + str(round(np.average(df_sorted.Contrast),3))

#legend
plt.legend(leg, shadow=True, fancybox=True, numpoints=1)






# differentiation 

fig2 = plt.figure(2)
fig2.set_size_inches(15.5, 8.5)

leg2 = ['']*(num_of_bins+1) # one theoretical plot + num_of_bins plots
list_of_dframes = np.array_split(halo_sorted, num_of_bins) 
for i in range (num_of_bins):
    
    df = list_of_dframes[i]
    halo_df_deriv = df.M200c
    nbins = 20
    # derivative dN/dlogM for each of the dframes
    dN = np.histogram(np.log10(halo_df_deriv), bins = nbins)[0]
    M = np.histogram(np.log10(halo_df_deriv), bins = nbins)[1]                  
    dlogM = np.diff(M)
    dN_dlogM = dN/dlogM
    #M_bin = [(10**M[i]+10**M[i+1])/2 for i in range (nbins)]
    M_bin =  [(M[i]+M[i+1])/2 for i in range (nbins)]
    M_bin =  [10**(M[i]) for i in range (nbins)]  
          
    log = 1
    if log==1:
       scale = 'logscale'
       plt.yscale('log')
       plt.xscale('log')
    else:
       scale = 'normal scale' 
    
    ttl = 'Mass function. Cube size='+str(maxcoord)+'Mpc. R='+str(round(R, 1))+'Mpc.'+' Scale='+scale
    plt.title(ttl)
       
    plt.xlabel('M, Msun/h')
    plt.ylabel('dN(>M)/dlogM')    
    
    plt.plot (M_bin, dN_dlogM/47104*num_of_bins, '.')
    
    leg2[i] = 'contr=' + str(round(np.average(df.Contrast),3))
    
# theory: 
M_theory = np.load('../../Press_Schechter/M_theory.npy')
dN_dlogM_theory = np.load('../../Press_Schechter/dN_dlogM_theory.npy')
plt.plot(M_theory, dN_dlogM_theory, '.')
leg2[num_of_bins] = 'Press-Schechter'
plt.legend(leg2, shadow=True, fancybox=True, numpoints=1)

plt.ylim(0, 1e8)

plt.show()
