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
R_str = '3.2'

# reading Volumes of bins
V_bins = np.load('./V_bins_nodes/V_bins_R_'+R_str+'_nodes.npy')
#V_bins = np.load('./V_bins/V_bins_R_'+R_str+'.npy')

#reading the halo DataFrame from the file
halo_df = pd.read_csv('./Halo_Dataframes/halo_df_mixed_R_'+R_str+'.csv', index_col=0)


halosize = np.size(halo_df, 0)

num_of_bins = 4

#dividing halos into bins
list_of_dframes = np.array_split(halo_df,   num_of_bins)   


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
       
    ttl = 'Number of halos with mass greater then M. '+str(num_of_bins)+' contrast bins. Cube size='+str(maxcoord)+'Mpc. R='+ R_str +'Mpc.'+' Scale='+scale
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
list_of_dframes = np.array_split(halo_df, num_of_bins) 

for i in range (num_of_bins):

    df = list_of_dframes[i]
    halo_df_deriv = df.M200c
    nbins = 10
    # derivative dN/dlogM for each of the dframes
    dN = np.histogram(np.log10(halo_df_deriv), bins = nbins)[0]
    M = np.histogram(np.log10(halo_df_deriv), bins = nbins)[1]                  
    dlogM = np.diff(M)
    dN_dlogM = dN/dlogM
    M_bin = [(10**M[s]+10**M[s+1])/2 for s in range (nbins)]
          
    log = 1
    if log==1:
       scale = 'logscale'
       plt.yscale('log')
       plt.xscale('log')
    else:
       scale = 'normal scale' 
    
    ttl = 'Mass function. Cube size='+str(maxcoord)+'Mpc. R='+R_str+'Mpc.'+' Scale='+scale
    plt.title(ttl)
       
    plt.xlabel('M, Msun/h')
    plt.ylabel('dN(>M)/dlogM')

    dN_dlogM = dN_dlogM/V_bins[i]
    plt.plot (M_bin, dN_dlogM)
    
    leg2[i] = 'contr=' + str(round(np.average(df.Contrast),3))
    
# theory: 
M_theory = np.load('../Press_Schechter/M_theory_cut.npy')
dN_dlogM_theory = np.load('../Press_Schechter/dN_dlogM_theory_cut.npy')
plt.plot(M_theory, dN_dlogM_theory)
leg2[num_of_bins] = 'Press-Schechter'
plt.legend(leg2, shadow=True, fancybox=True, numpoints=1)

plt.xlim(2.4e10, 3e14)
#plt.ylim(1e-4, 1e0)

plt.show()
