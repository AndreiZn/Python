#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:46:24 2017

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

def cll (c, h):
    #coord(i) is in space (m1,m2,m3)
    m1 = int(c[0]/h)
    m2 = int(c[1]/h)
    m3 = int(c[2]/h)       
    return ([m1,m2,m3])    

def p_in_cell (a1, a2, a3):
    index_array = []
    current_ind = linkedspace[a1][a2][a3]
    while (current_ind != -1):
        index_array.append(current_ind)
        current_ind = linkedl[current_ind]
    return index_array

t0 = time()
        
#number of cells and step
maxcoord = 32
numofcells = 20
hx, hy, hz = maxcoord/numofcells, maxcoord/numofcells, maxcoord/numofcells
h = hx #now hx = hy = hz = h

#max and min coordinates
xmin, ymin, zmin = 0, 0, 0
xmax, ymax, zmax = maxcoord, maxcoord, maxcoord

#getting coordinates
#coord = pgr.readsnap ('../rsnap_002.dat', 'pos', 1)
#nop = int(1e6) 
nop = int(np.size(coord)/3) #number of particles
#nop = int(512**3)
#coord = maxcoord*np.random.rand(nop,3)

#reading the halo list
prtclmass = 2.08108e07
halo_df = pd.read_csv('./32Mpc/halo_selected.csv', index_col = 0)
halo_df.sort_values(['M200c'], inplace=True)
halo_df = halo_df.loc[halo_df['M200c'] > 1e3*prtclmass]
int_array = np.array([i for i in range(np.size(halo_df, 0)) ])
halo_df.set_index(pd.Series(int_array, index = halo_df.index), inplace=True)
#random halo_df
halo_df.set_value(int_array, ['X', 'Y', 'Z'], maxcoord*np.random.rand(np.size(int_array),3));                 
halosize = int(np.size(halo_df)/4)
halo_arr = np.array(halo_df)


#values
R = 1.999*hx
V = 4/3*pi*R**3 
Mpc = 3.085678e24 #in sm
M_S = 1.989e33 #in g
UnitMass = 1e10*M_S 

#parameter
r1, r2, r3 = int(R/hx) + 1, int(R/hy) + 1, int(R/hz) + 1
r = r1 #now r1 = r2 = r3 = r

#creating the linkedspace
mx, my, mz = int((xmax - xmin)/h), int((ymax - ymin)/h), int((zmax - zmin)/h)
linkedspace = np.zeros((mx,my,mz), dtype=np.int) - 1

#to save coord[i] in a loop
part = np.zeros(3);

#coord = ([1.5*h, 1.5*h, 1.5*h], [1.5*h, 1.9*h, 1.5*h], [2.5*h, 1.9*h, 1.5*h], [1.5*h, 1.8*h, 1.5*h], [2.6*h, 1.9*h, 1.5*h])

#to save density of spheres encircling halos
halo_dens = np.zeros(halosize, dtype=np.int)

#index_arr will be used to save indexes of particles
index_arr = np.zeros(nop, dtype=np.int)

#cl will be used to save cell-coordinates
cl = np.array([0, 0, 0]);

t1 = time() - t0
print ('time to assign values to variables = ', t1)

#linked list
linkedl = np.zeros(nop, dtype=np.int)

ccode = \
        """
        int prtcl, m0, m1, m2, temp;
        for (prtcl = 0; prtcl<nop; prtcl = prtcl + 1)
        {
            cl(0) = coord(prtcl,0); cl(1) = coord(prtcl,1); cl(2) = coord(prtcl,2);
            m0 = int(cl(0)/h); m1 = int(cl(1)/h); m2 = int(cl(2)/h);
            temp = int(linkedspace(m0, m1, m2));
            linkedspace(m0,m1,m2) = prtcl;
            linkedl(prtcl) = temp;        
        }
        """
        
weave.inline (ccode, ['cl', 'nop', 'coord', 'linkedspace', 'linkedl', 'h'], type_converters = converters.blitz, compiler = 'gcc')  


t2 = time() - t0 - t1
print ('time to create and fill in the linked list = ', t2)

#calculating densities in spheres encircling halos

ccode2 = \
        """
        int halo_ind, ind, s, m, n, p, i, j, k, i_eff, j_eff, k_eff, current_ind;
        int counter = -1;
        double x, y, z;
        
        for(halo_ind = 0; halo_ind<halosize; halo_ind = halo_ind + 1)
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
                         
                         current_ind = linkedspace(i_eff, j_eff, k_eff);
                         counter = -1;
                         while (current_ind != -1){
                             counter+=1;
                             index_arr(counter) = current_ind;                             
                             current_ind = linkedl(current_ind);
                         }
                         
                         for (s = 0; s<=counter; s = s + 1){
                             ind = index_arr(s);
                             part(0) = coord(ind,0); part(1) = coord(ind,1); part(2) = coord(ind,2);
                             
                             if (i<0)           {part(0)-=maxcoord;}             
                             if (j<0)           {part(1)-=maxcoord;}         
                             if (k<0)           {part(2)-=maxcoord;}
                             if (i>=numofcells) {part(0)+=maxcoord;}
                             if (j>=numofcells) {part(1)+=maxcoord;}
                             if (k>=numofcells) {part(2)+=maxcoord;}
                             
                             if (((part(0) - x)*(part(0) - x) + (part(1) - y)*(part(1) - y) + (part(2) - z)*(part(2) - z)) <= R*R) {
                                 halo_dens(halo_ind)+=1;                             
                             }                                                         
                         }
                     } 
                 }
             }        
        }
        """
        
weave.inline (ccode2, ['halosize','halo_arr', 'halo_dens', 'coord', 'R', 'h', 'r', 'numofcells', 'linkedspace', 'linkedl', 'index_arr', 'part', 'maxcoord'], type_converters = converters.blitz, compiler = 'gcc')
        
halo_df.loc[:, 'Density'] = pd.Series(halo_dens, index=halo_df.index) 
t3 = time() - t0 - t1 - t2
print ('time to calculate density = ', t3)

#for halo_ind in range(halosize):
#    
#    cur_data = halo_df.iloc[halo_ind]
#    x, y, z = cur_data.X, cur_data.Y, cur_data.Z
#    m, n, p = int(x/h), int(y/h), int(z/h)
#    for i in range (m-r, m+r+1):
#        for j in range (n-r, n+r+1):
#            for k in range (p-r, p+r+1):
#                i_eff, j_eff, k_eff = i, j, k
#                if (i<0):           i_eff+=numofcells             
#                if (j<0):           j_eff+=numofcells             
#                if (k<0):           k_eff+=numofcells 
#                if (i>=numofcells): i_eff-=numofcells
#                if (j>=numofcells): j_eff-=numofcells
#                if (k>=numofcells): k_eff-=numofcells
#                ps = p_in_cell (i_eff,j_eff,k_eff)
#                
#                for ind in ps:
#                    part = list(coord[ind])
#                    
#                    if (i<0):           part[0]-=maxcoord             
#                    if (j<0):           part[1]-=maxcoord            
#                    if (k<0):           part[2]-=maxcoord 
#                    if (i>=numofcells): part[0]+=maxcoord
#                    if (j>=numofcells): part[1]+=maxcoord
#                    if (k>=numofcells): part[2]+=maxcoord
#                    
#                    if ((part[0] - x)**2 + (part[1] - y)**2 + (part[2] - z)**2 )<R**2:
#                        halo_df.loc[halo_ind, 'Density'] += 1