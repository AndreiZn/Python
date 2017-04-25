#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 03:47:09 2017

@author: andreyzn
"""

from __future__ import division #so that 3/2 = 1.5, but not 1 as in Python2.* 
import pygadgetreader as pgr
import numpy as np
from numpy import pi
from time import time
from scipy import weave
import pandas as pd
from scipy.weave import converters

t0 = time()
        
#number of cells and step
maxcoord = 64
numofcells = 40
#numofcells = 640
hx, hy, hz = maxcoord/numofcells, maxcoord/numofcells, maxcoord/numofcells
h = hx #now hx = hy = hz = h

#values
R = 1.999*hx
#R = 31.999*hx
V = 4/3*pi*R**3 
Mpc = 3.085678e24 #in sm
M_S = 1.989e33 #in g
UnitMass = 1e10*M_S 

#max and min coordinates
xmin, ymin, zmin = 0, 0, 0
xmax, ymax, zmax = maxcoord, maxcoord, maxcoord

#reading the halo list
prtclmass = 2.08108e07
halo_df = pd.read_csv('./Halo dataframes for different radii/halo_df_64Mpc_R_3.2.csv', index_col = 0)

#normalized density:
temp_arr = np.array(halo_df.Density)/V
halo_df.loc[:, 'Density'] = pd.Series(temp_arr, index=halo_df.index)

#Contrast column
aver_dens = 512**3/maxcoord**3
halo_df.loc[:, 'Contrast'] = (halo_df.Density - aver_dens)/aver_dens

#parameter
r1, r2, r3 = int(R/hx) + 1, int(R/hy) + 1, int(R/hz) + 1
r = r1 #now r1 = r2 = r3 = r

t1 = time() - t0
print ('time to assign values to variables = ', t1)

num_of_bins = 4

#sorting the halo DataFrame
halo_sorted = halo_df.sort_values(['Contrast'])

#Index is mixed up. Changing it to 0, 1, 2, ... int_array
int_array = np.array([i for i in range(np.size(halo_sorted, 0)) ])
halo_sorted.loc[:,'new_index'] = pd.Series(int_array, index = halo_sorted.index)
halo_sorted.set_index('new_index', inplace=True)

#dividing halos into bins
list_of_dframes = np.array_split(halo_sorted, num_of_bins)   

#array with bins' volumes
bins_V = np.array([0]*num_of_bins, dtype=float)
cell_V = h**3

for bin_ind in range (num_of_bins):

    df = list_of_dframes[bin_ind]
    halo_arr = np.array(df)
    halos_in_bin = int(np.size(df)/6) #6 columns
    #creating the linkedspace
    mx, my, mz = int((xmax - xmin)/h), int((ymax - ymin)/h), int((zmax - zmin)/h)
    space = np.zeros((mx,my,mz), dtype=np.int)      
            
    #find cells covered by spheres with centers at (x,y,z), where (x,y,z) are coordinates of halos
    
    ccode = \
        """
        int halo_ind, ind, m, n, p, i, j, k, i_eff, j_eff, k_eff;
        double x, y, z;
        
        for(halo_ind = 0; halo_ind<halos_in_bin; halo_ind = halo_ind + 1)
        {
             x = halo_arr(halo_ind, 0); y = halo_arr(halo_ind, 1); z = halo_arr(halo_ind, 2);
             m = int(x/h); n = int(y/h); p = int(z/h);
             for (i = m - r; i<m+r+1; i = i + 1) {
                 for (j = n - r; j<n+r+1; j = j + 1){
                     for (k = p - r; k<p+r+1; k = k + 1){ 
                         i_eff = i; j_eff = j; k_eff = k;
                         if (i<0)           {i_eff+=numofcells;}             
                         if (j<0)           {j_eff+=numofcells;}             
                         if (k<0)           {k_eff+=numofcells;} 
                         if (i>=numofcells) {i_eff-=numofcells;}
                         if (j>=numofcells) {j_eff-=numofcells;}
                         if (k>=numofcells) {k_eff-=numofcells;}
                         
                         if (space(i_eff, j_eff, k_eff) == 0) {                              
                                 space(i_eff, j_eff, k_eff) = 1;
                                 bins_V(bin_ind)+=cell_V;
                                 
                         }
                                                                                      
                         
                     } 
                 }
             }        
        }
        """                                               
                                                   
    weave.inline (ccode, ['halos_in_bin','halo_arr', 'R', 'h', 'r', 'numofcells', 'space', 'maxcoord', 'bin_ind', 'bins_V', 'cell_V'], type_converters = converters.blitz, compiler = 'gcc')
    
    
t2 = time() - t0 - t1
print ('time to calculate volumes of bins = ', t2)   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    