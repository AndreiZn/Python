# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 17:16:18 2016

@author: Andrey
"""
from __future__ import division #so that 3/2 = 1.5, but not 1 as in Python2.* 
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi as pi
import scipy.integrate as integrate

A = 4.375e6
Omega = 0.27
h = 0.73
M_S = 1.98e30
Mpc = 3.0856776e22
#rho = 0.3*9.31e-27*Mpc**3/h**2
rho = 1.88e-26*Mpc**3 #in h**2 units 
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
    return ( integrate.quad(lambda t: np.exp(3*t)*P(t)*(W(t,R))**2/(2*pi**2), np.log(1e-2), np.log(1e5))[0] )  

def der_sgm2(M):         
    return ( integrate.quad(lambda t: np.exp(3*t)*P(t)*2*br(t,M)*(br(t,M)-al(t,M)**2*np.sin(al(t,M)))/(al(t,M)**6*M*2*pi**2), np.log(1e-2), np.log(1e5))[0] )

def dersgm(M):
    return der_sgm2(M)/(2*np.sqrt(sgm2(M)))

def n(M):
    return np.sqrt(2/pi)*rho*dersgm(M)*sgm_s*np.exp(-sgm_s**2/(2*sgm2(M)))/(M*sgm2(M))
    
M_0 = M_S*1e5
M_F = M_S*1e15
M = np.logspace(np.log10(M_0), np.log10(M_F), num = 200)
simga2 = [sgm2(M[i])/h**2 for i in range (len(M))]
dN_dM = [n(M[i]) for i in range (len(M))]
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
M_x = M*h/M_S
plt.plot (M_x, dN_dlogM, '.')

plt.show()

# save plotted arrays
np.save('M_theory', M_x)
np.save('dN_dlogM_theory', dN_dlogM)