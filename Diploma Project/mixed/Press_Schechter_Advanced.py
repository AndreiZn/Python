# -*- coding: utf-8 -*-
"""
Created on Thu May 18 00:50:38 2017

@author: Andrey

"""
from __future__ import division #so that 3/2 = 1.5, but not 1 as in Python2.* 
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi as pi
import scipy.integrate as integrate
import pandas as pd

A = 4.375e6
#Omega = 0.27
Omega = 0.307115
#h = 0.73
h = 0.6777
M_S = 1.98e30 # in kg
Mpc = 3.0856776e22 # in meters
#rho = 0.3*9.31e-27*Mpc**3/h**2
rho = 1.88e-26*Mpc**3/M_S #in h**2 units 
sgm_s=1.686
#k_0 = Omega*h**3/(3.086e22)
k_0 = Omega*h**2
g = 4*pi*rho/3

def Rad (M):
    return (M/g)**(1/3)     
    
def P(x):
    k = np.exp(x)
    return A*k*(T(k))**2
            
def T(x):
    q = x/(Omega*h**2)
    return np.log(1+2.34*q)*(2.34*q)**(-1)*(1+3.89*q+(16.1*q)**2+(5.64*q)**3+(6.71*q)**4)**(-1/4)

def W(x,R):
    k = np.exp(x)
    return (3*np.sin(k*R)-3*k*R*np.cos(k*R))/(k*R)**3
    #return 1/(R)**4/k**3    

def al(t, M):
    return np.exp(t)*g**(-1/3)*M**(1/3)
    
def br(t, M):
    a = al(t, M)
    return 3*np.sin(a) - 3*a*np.cos(a)
    
def sgm2(M):
    R = Rad (M)
    #return ( integrate.quad(lambda x, a: x**2*P(x)*(W(x,a))**2/(2*pi**2) , 0 , 1e3, args=(R,))[0] )
    #return ( integrate.quad(lambda t: np.exp(3*t)*P(t)*(W(t,R))**2/(2*pi**2), np.log(1e-3*h/Mpc), np.log(1e1*h/Mpc))[0] )    
    return ( integrate.quad(lambda t: np.exp(3*t)*P(t)*(W(t,R))**2/(2*pi**2), np.log(1e-5), np.log(1e5))[0] )  

def der_sgm2(M):         
    return ( integrate.quad(lambda t: np.exp(3*t)*P(t)*2*br(t,M)*(br(t,M)-al(t,M)**2*np.sin(al(t,M)))/(al(t,M)**6*M*2*pi**2), np.log(1e-2), np.log(1e5))[0] )

def dersgm(M):
    return der_sgm2(M)/(2*np.sqrt(sgm2(M)))

def n(M):
    return np.sqrt(2/pi)*rho*dersgm(M)*sgm_s*np.exp(-sgm_s**2/(2*sgm2(M)))/(M*sgm2(M))
    

#R_str_array = ['1.0', '2.0', '2.6', '3.2', '4.0', '6.4', '8.0']
#R_str_array = ['4.0', '6.4', '8.0']
R_str_array = ['4.0']
for k in range (np.size(R_str_array)):
    
    R_str = R_str_array[k]
    
    halo_df = pd.read_csv('../Halo_Programs/Halo_Dataframes/halo_df_mixed_cut_R_'+R_str+'.csv', index_col=0)
    
    num_of_bins = 4
    
    #dividing halos into bins
    list_of_dframes = np.array_split(halo_df,   num_of_bins)   
    
    for p in range (num_of_bins):
        df = list_of_dframes[p]
        delta_min = np.min(df.Contrast)
        sgm_s = 1.686 - delta_min
        #delta_max = np.max(df.Contrast)
        #sgm_s = 1.686 + delta_max
        M_0 = 1e10*h
        M_F = 1e15*h
        M = np.logspace(np.log10(M_0), np.log10(M_F), num = 200)
        simga2 = [sgm2(M[i]) for i in range (len(M))]
        
        #if sgm_s > 0:
        dN_dM = [n(M[i]) for i in range (len(M))]
        #else:
        #    dN_dM = [-n(M[i]) for i in range (len(M))]
            
        dN_dlogM = [dN_dM[i]*M[i] for i in range (len(M))]
        
        fig2 = plt.figure(2)
        
        mode = 'MFunc'
        log = 1
        if mode=='MFunc':
            plt.title('The Press-Schechter mass function')
            plt.xlabel('M, M_Sun/h')
            plt.ylabel('dN(>M)/dlogM')
        if log == 1:    
            plt.xscale('log')
            plt.yscale('log')
        
        #print('!')        
        plt.plot (M, dN_dlogM, '.')
        
        # save plotted arrays
        np.save('npy_arrays_4/'+R_str+'/M_theory_adv_'+R_str+'_'+str(round(sgm_s, 1)), M)
        np.save('npy_arrays_4/'+R_str+'/dN_dlogM_theory_adv_'+R_str+'_'+str(round(sgm_s, 1)), dN_dlogM)
        
    plt.show()    