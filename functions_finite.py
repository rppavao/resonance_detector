# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 13:19:16 2018

@author: rpavao
"""
from cmath import phase,sqrt, exp,log, pi
#import numpy as np
#import numpy.linalg as nplin
#import scipy as sc
#import scipy.integrate as scint
#import time
#import matplotlib.pyplot as plt
#import itertools as itr
#from mpl_toolkits.mplot3d import Axes3D
#from numba import jit
from matrix import masses
import config

def Ef(M,m,z):
    return (z - m**2. + M**2.)/(2.*sqrt(z))


#


#Energy of Baryon in function of s

#def Vf1(z,x,y):
#	return (sqrt((sqrt(z) + MmB[x-1])/MmB[x-1])*(2.*sqrt(z) - MmB[x-1] - MmB[y-1])* sqrt((sqrt(z) + MmB[y-1])/MmB[y-1])*Dm[x-1]#[y-1])/(8.*Mf[x-1]*Mf[y-1])


def sf(m,M,nc):
	if nc=='+':
		return (m+M)**2.
	if nc=='-':
		return (m-M)**2.
	else:
		return 'ERROR'

def rhof(z,m,M,nc):
	if nc=='+':
		return abs(z-sf(m,M,'+'))
	if nc=='-':
		return abs(z-sf(m,M,'-'))
	else:
		return 'ERROR'	

def argf(z,nc):
	x=phase(z)
	if nc == 0:
		if x < 0:
			return x+2.*pi
		else:
			return x
	if nc== 1:
		if x==pi:
			return -pi
		else:
			return x
	else:
		return 'ERROR'


def Rf(z,n,m,M):
	return (sqrt(rhof(z, m, M, '+'))*exp (1.j*argf(z - sf (m, M, '+'), 0)/2.) + sqrt(rhof(z, m, M, '-'))* \
 exp(1.j*argf(z - sf (m, M, '-'), 1.)/2.)*exp(1.j*n*pi))/(sqrt(rhof(z, m, M, '+'))*exp(1.j*argf (z - sf(m, M,'+'), 0)/2.) - \
 sqrt(rhof(z, m, M, '-'))*exp(1.j*argf(z - sf(m, M, '-'), 1.)/2.)*exp(1.j*n*pi))

def Lf(z, n, m, M):
	return (sqrt(rhof(z, m, M, '+')*rhof (z, m, M, '-'))/z) * exp(1.j*(argf(z - sf (m, M, '+'), 0) + \
  argf(z - sf (m, M,'-'), 1) + 2.*n*pi)/2.)*(log(abs (Rf (z, n, m, M))) + 1.j*argf(Rf(z, n, m, M), 0) - 2.*pi*1.j)
#============================================================================== 
"""
Important functions that the code must always have:
"""
class Error_Index(Exception):
    pass
##Define the Loop function
def G(z, i,*par,LnSRS,box=False,l):
    class Error_In_Box(Exception):
        pass 
    
    try:
        while len(par)==1:
            par=par[0]
    except TypeError:
        pass
#    par = par[0][0][0]
    if i-1 <0:
        raise Error_Index("The indices start at 1.")        
    xdc=par[0]
    cut=par[1]
    J=par[3]
    Mlit=masses(par[2],2,0.5)
    Mth2=Mlit[0]**2.#+Mlit[1]**2.
    alfa = 1.
    mu = sqrt(alfa*Mth2)
    
    Mass_tuple=masses(par[2],i-1,J)

    M=Mass_tuple[0].real
    m=Mass_tuple[1].real
#    LnSRS=nSRS(J)
    	#Vector n to specify the riemann sheets

    if ((sqrt(z)).real)**2. < (m+M)**2.: ################
    	n=LnSRS[0][i-1]
    else:
        if z.imag>0:
    	    n=LnSRS[0][i-1]
        else:
            n=LnSRS[1][i-1]
        
    if box is False:           
        Gdim=(2.*M/(4.*pi)**2.)*(((M**2. - m**2.)/mu**2. - (M - m)/(M + m))*log(float(M/m)) + Lf(mu**2., 0, m, M))
        Gcut=M*(m*log(m/(cut+sqrt(cut**2+m**2)))+M*log(M/(cut+sqrt(cut**2+M**2))))/(4.*pi**2*(m+M))
        
        return (2.*M/(4.*pi)**2.)*(((M**2. - m**2.)/z - (M - m)/(M + m))*log(float(M/m)) + Lf(z, n, m, M))-\
            (1-xdc)*Gdim+xdc*Gcut
    elif box is True:    
        import os,sys
        import Gfinite

        def Gbox(e,m1,m2,mu,l):
              lfv=    l
              return Gfinite.gfvol(e,m1,m2,lfv,mu)
          
        Gdim=(2.*M/(4.*pi)**2.)*(((M**2. - m**2.)/mu**2. - (M - m)/(M + m))*log(float(M/m)) + Lf(mu**2., 0, m, M))
        Gcut=M*(m*log(m/(cut+sqrt(cut**2+m**2)))+M*log(M/(cut+sqrt(cut**2+M**2))))/(4.*pi**2*(m+M))        
        Gloop=(2.*M/(4.*pi)**2.)*(((M**2. - m**2.)/z - (M - m)/(M + m))*log(float(M/m)) + Lf(z, 0, m, M))-\
            (1-xdc)*Gdim+xdc*Gcut            
            
        return Gbox(sqrt(z),m,M,mu,l)+Gloop.real
    else:
        raise Error_In_Box("Use only True or False for box parameter")

'''
G function for finite space
'''
      

#Element x,y of the matrix V of interaction
def Vf(z,x,y,*par):
    '''
    This gives V_xy, x and y start at 1! 
    '''
    try:
        while len(par)==1:
            par=par[0]
    except TypeError:
        pass
    
    if x-1 <0 or y-1 <0:
        raise Error_Index("The indices start at 1.")
    J=par[3]
    Mass_tuplex=masses(par[2],x-1,J)
    Mass_tupley=masses(par[2],y-1,J)
    
    Mi=Mass_tuplex[0]
    mi=Mass_tuplex[1]
    Mo=Mass_tupley[0]
    mo=Mass_tupley[1]
    Mfx=Mass_tuplex[2]
    Mfy=Mass_tupley[2]
    Dm=Mass_tuplex[3]
#    print(Dm)
#    print(Mi)
#    print(mi)
#    print(x,y)
#    print((2.*sqrt(z) - Mi - Mo)*Dm[x-1][y-1]/(4.*Mfx*Mfy))
#    stop
#    config.message.add("Relativistic Corrections")
#    return (sqrt((Ef(Mi,mi, z) + Mi)/Mi)*(2.*sqrt(z) - Mi - Mo)* \
#            sqrt((Ef(Mo, mo, z) + Mo)/Mo)*Dm[x-1][y-1])/(8.*Mfx*Mfy)
    return (2.*sqrt(z) - Mi - Mo)*Dm[x-1][y-1]/(4.*Mfx*Mfy)


