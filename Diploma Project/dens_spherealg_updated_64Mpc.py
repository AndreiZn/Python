#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 00:43:18 2017

@author: andreyzn

This version of the program has the main body completely rewritten in C
"""



from __future__ import division #so that 3/2 = 1.5, but not 1 as in Python2.* 
import pygadgetreader as pgr
import numpy as np
from numpy import pi
from time import time
from scipy import weave
from scipy.weave import converters

t0 = time()
coord = pgr.readsnap ('snapshot_009', 'pos', 1) 
t1 = time() - t0
maxcoord = 64
numofcells = 20
hx, hy, hz = maxcoord/numofcells, maxcoord/numofcells, maxcoord/numofcells
h = hx #now hx = hy = hz = h
numofsph = numofcells + 1 # = number of knots
R = 1.999*hx
V = 4/3*pi*R**3 
Mpc = 3.085678e24 #in sm
M_S = 1.989e33 #in g
UnitMass = 1e10*M_S 

#nop = int(np.size(coord)/3) #number of particles
nop = int(1e7);

xmax, ymax, zmax = 32, 32, 32
xmin, ymin, zmin = 0, 0, 0
part = np.zeros(3)
density = np.zeros((numofsph,numofsph,numofsph))

r1 = int(R/hx) + 1
r2 = int(R/hy) + 1
r3 = int(R/hz) + 1
r = r1 #now r1 = r2 = r3 = r

def boundcond ():
    lstsph = numofsph - 1
    density[lstsph,:,:] = density [0, :, :]
    density[:,lstsph,:] = density [:, 0, :]
    density[:,:,lstsph] = density [:, :, 0]

ccode = \
        """
        int prtcl = 0;
        int i, j, k, m, n, p, meff, neff, peff;
        double xsp, ysp, zsp;
        for (prtcl = 0; prtcl<nop; prtcl = prtcl +1)
        {
            part(0) = coord(prtcl,0); part(1) = coord(prtcl,1); part(2) = coord(prtcl,2);
            m = int(part(0)/h); n = int(part(1)/h); p = int(part(2)/h);
        
            for (i = -r + 1; i < r + 1; i = i + 1){
                for (j = -r + 1; j < r + 1; j = j + 1){
                    for (k = -r + 1; k < r + 1; k = k + 1){
                        xsp = (m+i)*h; ysp = (n+j)*h; zsp = (p+k)*h; 
                        if (((part(0)-xsp)*(part(0)-xsp) + (part(1)-ysp)*(part(1)-ysp) + (part(2)-zsp)*(part(2)-zsp)) <= R*R) {
                            meff = m + i; neff = n + j; peff = p + k;
                            if (meff<0) {meff+=numofsph-1;}
                            if (neff<0) {neff+=numofsph-1;}
                            if (peff<0) {peff+=numofsph-1;}
                            if (meff>=numofsph) {meff-=numofsph-1;}
                            if (neff>=numofsph) {neff-=numofsph-1;}
                            if (peff>=numofsph) {peff-=numofsph-1;}
                            density (meff, neff, peff)+=1;
                        }
        }
        }
        }
        }
        """
weave.inline (ccode, ['density', 'nop', 'numofsph', 'part', 'R', 'h', 'r', 'coord'], type_converters = converters.blitz, compiler = 'gcc')   
density = density/V;                        
boundcond()

t2 = time() - t1 - t0
print ('time for reading = ', t1, 'time for calculating density', t2)  
 
#index = peff+neff*numofsph+meff*numofsph*numofsph;
                            #density[index]+=1;               