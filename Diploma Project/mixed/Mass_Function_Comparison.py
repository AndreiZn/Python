# -*- coding: utf-8 -*-
"""
Created on Thu May 18 01:34:26 2017

@author: Andrey
"""

from __future__ import division #so that 3/2 = 1.5, but not 1 as in Python2.* 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab as plb

maxcoord = 64

R_str_array = ['1.0', '2.0', '2.6', '3.2', '4.0', '6.4', '8.0']
#R_str_array = ['4.0', '6.4', '8.0']
#R_str_array = ['3.2', '8.0']
#R_str_array = ['4.0']

for k in range (np.size(R_str_array)):
    
    R_str = R_str_array[k]
    
    # reading Volumes of bins
    #V_bins = np.load('./V_bins_nodes/V_bins_R_'+R_str+'_nodes.npy')
    V_bins = np.load('./V_bins/V_bins_R_'+R_str+'.npy')
    #V_bins = V_bins
    
    #reading the halo DataFrame from the file
    halo_df = pd.read_csv('./Halo_Dataframes/halo_df_mixed_cut_R_'+R_str+'.csv', index_col=0)
    
    
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
    
    fig2 = plt.figure(k+2)
    fig2.set_size_inches(15.5, 8.5)
    
    leg2 = ['']*(num_of_bins*2) # one theoretical plot + num_of_bins plots
    list_of_dframes = np.array_split(halo_df, num_of_bins) 
    
    for i in range (num_of_bins):
    
        df = list_of_dframes[i]
        #delta_min = np.min(df.Contrast)
        #sgm_s = 1.686 + delta_min
        delta_min = np.min(df.Contrast)
        sgm_s = 1.686 - delta_min
        #delta_max = np.max(df.Contrast)
        #sgm_s = 1.686 + delta_max
        
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
        
        leg2[2*i] = 'Simulation, Contrast = ' + str(round(np.min(df.Contrast),3)) + ' : ' + str(round(np.max(df.Contrast),3))
        
        # theory: 
        M_theory = np.load('../Press_Schechter/npy_arrays_min2/'+R_str+'/M_theory_adv_'+R_str+'_'+str(round(sgm_s, 1))+'.npy')
        dN_dlogM_theory = np.load('../Press_Schechter/npy_arrays_min2/'+R_str+'/dN_dlogM_theory_adv_'+R_str+'_'+str(round(sgm_s, 1))+'.npy')    
        #M_theory = np.load('../Press_Schechter/npy_arrays_3/M_theory_adv_'+R_str+'_'+str(round(sgm_s, 1))+'.npy')
        #dN_dlogM_theory = np.load('../Press_Schechter/npy_arrays_3/dN_dlogM_theory_adv_'+R_str+'_'+str(round(sgm_s, 1))+'.npy')
        plt.plot(M_theory, dN_dlogM_theory, '.')
        leg2[2*i+1] = 'Press-Schechter, contr = ' + str(round(delta_min,3)) + ' : ' + str(round(np.max(df.Contrast),3))
        #leg2[2*i+1] = 'Press-Schechter, delta_c + ' + str(round(delta_max,1))
    
    plt.legend(leg2, shadow=True, fancybox=True, numpoints=1, prop={'size':9})
    
    plt.xlim(2.4e10, 1e13)
    plt.ylim(1e-4, np.max(dN_dlogM_theory))
    
    plt.show()
    
    plb.savefig('MF_Comparison_'+R_str+'.png', bbox_inches='tight')
    