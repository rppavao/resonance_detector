#from cmath import phase,sqrt, exp,log, pi
import numpy as np
import numpy.linalg as nplin
#import scipy as sc
#import scipy.integrate as scint
#import time
#import matplotlib.pyplot as plt
import itertools as itr
#from mpl_toolkits.mplot3d import Axes3D
#from numba import jit
#import matrix
from matrix import dim
#import functions
from functions import Vf, G
import config
from importlib import reload
#Here we define the lightest channel
def Tmat(z,*par,LnSRS,box=False):
    try:
        while len(par)==1:
            par=par[0]
    except TypeError:
        pass
    J=par[3]
    cdim=dim(J)
    Vm=[[Vf(z, i, j,par) for j in np.arange(1,cdim+1)] for i in np.arange(1,cdim+1)]

    Gm=[[ 0  for j in np.arange(1,cdim+1)] for i in np.arange(1,cdim+1)]
    for i in np.arange(cdim):
        Gm[i][i]=G(z,i+1,par,LnSRS=LnSRS,box=box)

    Tm1=nplin.inv(Vm)-Gm
    Tm=nplin.inv(Tm1)
      
#    Tm=max([sum(Tm2[i]) for i in np.arange(len(Tm2[0]))])
    return Tm

def T2(z,channel1=None,channel2=None,*par,LnSRS,box=False):
    if channel1 is not None and channel2 is not None:
        config.message.add('NOTE: We are using T['+str(channel1)+']['+str(channel2)+']')
        Tm3=abs(Tmat(z,par,LnSRS=LnSRS,box=box)[channel1-1][channel2-1])
    else:
        Tm2=abs(Tmat(z,par,LnSRS=LnSRS,box=box))
        Tm3=max([sum(Tm2[i]) for i in np.arange(len(Tm2[0]))])
    return Tm3

def getlists(Emin,Emax,Wmin,Wmax,Ep,Wp):
    counter = itr.count()
    X=[i*float((Emax-Emin)*1.0/Ep)+Emin for i in np.arange(Ep)]
    Y=[j*float((Wmax-Wmin)*1.0/Wp)+Wmin for j in np.arange (Wp)]
    return (X,Y)
#       
#List0=jit(getlists)
def getpoles(Emin,Emax,Wmin,Wmax,Ep,Wp,ch1,ch2,*par,LnSRS,box=False):
    print(par)
    List1=getlists(Emin,Emax,Wmin,Wmax,Ep,Wp)
    X=List1[0]
    Y=List1[1]
    XS,YS = np.meshgrid(X,Y)
    maxlist=[]
    Z=np.array([[abs(T2((XS[i][j]-YS[i][j]*0.5j)**2,ch1,ch2,par,LnSRS=LnSRS,box=box)) for j in range(len(XS[i]))]for i in range(len(XS))])
    for i in Z:
         maxlist.append(max(i))
    zmax=max(maxlist)
    return (XS,YS,Z/zmax,zmax)

