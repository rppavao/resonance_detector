#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:00:40 2019

@author: rpavao
"""
import time as t
import numpy as np
import numpy.linalg as nplin
from polescalc import Tmat
from functions import G,Vf
from cmath import sqrt,pi
from scipy import optimize as op
import matplotlib.pyplot as plt
from plot_functions import get_xyz_function,change_scale,get_xyz_files
from matrix import MmB,Mmm,dim
from poles import HiddenPrints
import config
from numba import jit
if config.message != set():
    print(config.message)
    
def Efinite(m1,m2,l,n):
    q2=(2*pi/l)**2*n
    return sqrt(q2+m1**2)+sqrt(q2+m2**2)    

'''
Parameters
'''
direct= "./finite_data/"
name0="finite_ymass="
ti= t.time()
for ymass in range(1):
    parameters=(0,630,ymass,0.5)
    mpar=parameters[2]
    decimal=7
    fermi=197.326
    lin=1.5/fermi
    lfin=10./fermi
    #dl=.5/fermi
    nlpoints=30
    nmax=4
    npoles=10
    dEmax=1.
    Emin=min([Efinite(MmB[i][ymass],Mmm[i][ymass],lin,0).real for i in range(len(MmB))])-180
    Emax=Emin+550
    eps=10**(-6)
    #name="finite_"+str(lin)+"-"+str(lfin)
    name=direct+name0+str(ymass+1)   
    make_data=True
        
        
        
        
    '''    
    Define functions 
    '''
#    thresh=[round(Mmm[i][0]+MmB[i][0],decimal) for i in range(len(MmB))]
    
    
    RS=[[0 for i in range(dim(0.5))],[1 for i in range(dim(0.5))]]
    
    def delt(i,j):
        if i==j:
            return 1
        else:
            return 0
    
    def detT(E,l):
    #    l=l
        J=parameters[3]
        cdim=dim(J)
        Vm=[[Vf(E**2, i, j,parameters) for j in np.arange(1,cdim+1)] for i in np.arange(1,cdim+1)]
        
        Gm=[[ 0  for j in np.arange(1,cdim+1)] for i in np.arange(1,cdim+1)]
        for i in np.arange(cdim):
            Gm[i][i]=G(E**2,i+1,parameters,LnSRS=RS,box=True,l=l)
        Ident=np.array([[delt(i,j) for j in np.arange(1,cdim+1)] for i in np.arange(1,cdim+1)])
    #    return nplin.det(nplin.inv(Tmat(E**2,parameters,LnSRS=RS,box=True))).real
        return  nplin.det(Ident-np.matmul(np.array(Vm),np.array(Gm)))

    
    '''try to cify'''
    #@jit(nopython=True)
    #def findguess(XY):
    #    X=XY[0]
    #    Y=XY[1]c
    #    R=[]
    #    for i in range(1,len(Y)):
    #        if np.sign(Y[i-1]) != np.sign(Y[i]):
    #            R.append((X[i-1],X[i]))
    #    return R
    #
    #def roots(f,a,b,npoints=20):
    #    XY=get_xyz_function(f,a,b,npoints)
    #    R=findguess(XY)
    #    POLES=[]
    #    for i in R:
    #        Emin=op.bisect(f,i[0],i[1])
    #        if round(Emin,4) not in thresh:
    #            POLES.append(Emin)
    #    return POLES
                  
    def roots(f,a,b,dEmax):
    #    if dEmax>(b-a)/5:
    #        dEmax=(b-a)*0.01
        npoints=int(round((b-a)/dEmax,0))
    
        if dEmax>(b-a)/5:
            npoints=int(200)    
        
        
        XY=get_xyz_function(f,a,b,npoints)
        X=XY[0]
        Y=XY[1]
    #    print("E_KN_0=",Efinite(MmB[0][0],Mmm[0][0],6.5/fermi,0))
    #    print("E_Sp_1=",Efinite(MmB[1][0],Mmm[1][0],6.5/fermi,1))    
    #    print("a=",X[0],"\t b=",X[-1],"\t b-1=",X[-2])
    #    print("sign(a)=",Y[0],"\t sign(b)=",Y[-1],"\t sign(b-1)=",Y[-2],"\n")
        R=[]
        for i in range(1,len(Y)):
            if np.sign(Y[i-1]) != np.sign(Y[i]):
    #            print("X=",X[i-1].real,"\t",X[i].real)
    #            print("Y=",round(Y[i-1].real,4),"\t",round(Y[i].real,4))
                Emin=op.bisect(f,X[i-1],X[i])
    #            if round(Emin,decimal) not in thresh:
                R.append(Emin)
        return R
    
    
    
    def froot(l,E,X,npoles,dEmax,eps):
    #    m=totnpoints/(Emax-Emin)
        i=X.index(l*fermi)
        Enew=[]
        for c in range(len(E)):
            for n in range(len(E[c])):
                En=E[c][n][i]
                if En>= Emin  and En<= Emax:
                    Enew.append(E[c][n][i])
        Enew.sort()
        Enew.insert(0,Emin)
        Enew.append(Emax)
        R=[]
    #    print("LISTN=",Enew[:10])
        for e in range(1,len(Enew)):
    #        npoints=int(Enew[e+1]-Enew[e])
    #        if npoints==0:
    #            continue
    #        if Enew[e+1]-Enew[e] <1:
    #            npoints=
    #        print("a#eps=",Enew[e-1]+eps,"b#eps=",Enew[e]-eps)
    #        print("a#=",Enew[e-1],"b#=",Enew[e])
            R+=roots(lambda x: detT(x,l),Enew[e-1]+eps,Enew[e]-eps,dEmax)
    #        print("dEmax=",dEmax)
            if len(R)>npoles:
                break
        R.sort()
        while len(R) < npoles:
            R.append(10**30)
    #    return roots(lambda x: detT(x,l),Emin,Emax,npoles=npoles)
        return R 
    
    
    plt.figure(ymass+1)    
    Ethfin=[[[] for j in range(nmax+1)] for i in range(dim())]
    
    for n in range(nmax+1):
        for i in range(len(MmB)):    
            m1=MmB[i][mpar]
            m2=Mmm[i][mpar]
        #    print(m1+m2)
    #        print(lin,lfin,int(round((lfin+dl-lin)/dl,0)))
    #        nlpoints=int(round((lfin+dl-lin)/dl,0))
            XYZ=get_xyz_function(lambda l: Efinite(m1,m2,l,n).real,lin,lfin,nlpoints)
            Ethfin[i][n]=XYZ[1]
            Xfin=[x*fermi for x in XYZ[0]]
            plt.plot(Xfin,XYZ[1],color='r',linewidth=1.,linestyle='dashed')
    
    '''
    End of functions
    '''
    t1= t.time()    
    ##Emin,Emax=140,1750        
    #l=6.5/fermi
    #axes = plt.gca()
    #axes.set_ylim([-1,1])
    #axes.set_xlim([1425,1450])
    ##PONT=roots(lambda x: detT(x,l),1425,1450,0.1)
    #PONT2=froot(l,Ethfin,Xfin,5,1.,eps=10**(-6))
    #plt.axhline(y=0,color='k',linewidth=1.)
    #XY=get_xyz_function(lambda xvar: detT(xvar,l),Emin,Emax,10000)
    #X=XY[0]
    #Y=XY[1]
    #
    #i=Xfin.index(l*fermi)
    #Enew=[]
    #for c in range(len(Ethfin)):
    #    for n in range(len(Ethfin[c])):
    #        Enew.append(Ethfin[c][n][i])
    #Enew.sort()
    #
    ##import sys
    ##sys.exit("Error message")
    #plt.plot(X,Y)
    #for i in Enew[:7]:
    #    plt.axvline(x=i,color='k',linestyle='--')
    #for i in PONT2:
    #    plt.axvline(x=i,color='r',linestyle='-.')
    #plt.savefig("pedro.pdf")
    ##print(PONT)
    #t2= t.time()
    #
    #print("\nThe program took "+str(round(t2-t1,4))+" seconds to run.")
    #
    #import sys
    #sys.exit("Error message")
            
    if make_data==True:
        f=open(name,"w")
        l=lin
        dl=(lfin-lin)/nlpoints
        for i in range(nlpoints+1):
            l=lin+(lfin-lin)*i/nlpoints
            Ystr=''
            Y=froot(l,Ethfin,Xfin,npoles,dEmax,eps)
            for j in Y:
                Ystr+=str(j)+"\t"
    #        X=[l*fermi for i in Y]
        #    print(Ystr,Y)
            f.write(str(l*fermi)+"\t"+Ystr+"\n")
            print("l=",round(l*fermi,2))
    #        l+=dl
        #    plt.scatter(X,Y,color='k')  
        f.close()
        
    axes = plt.gca()
    axes.set_xlim([lin*fermi,lfin*fermi])
    axes.set_ylim([Emin,Emax-50])
    XYZ=get_xyz_files(name)
    X=XYZ[0]
    for Y in XYZ[1:]:
        plt.plot(X,Y,color='k',linewidth=1.)
    
    
    plt.savefig(name+".pdf")
    
    t2= t.time()
    
    print("\nThe program took "+str(round(t2-t1,4))+" seconds to run.")
    
    
tf= t.time()

print("\nThe program took a total of "+str(round((tf-ti)/60,4))+" minutes to run.")
