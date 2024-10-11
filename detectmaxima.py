# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 11:19:52 2018

@author: rpavao
"""

import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes as ax
import time as t

import sklearn 
from sklearn import cluster
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans

import importlib as ilib
import polescalc as pcalc
import config
from matrix import dim
from poles import pole, HiddenPrints

def cut(Z,zmax=None):
    if zmax is not None:
        new_Z=[]
        for i in Z:
            new_Z.append([])
            for j in i:
                if j>zmax:
                    newz=zmax
                else:
                    newz=j
                new_Z[-1].append(newz)
        return np.array(new_Z)
    else:
        return Z
    
def plotmax(file,count_del=0,masstf=False,graphcut=100):
    data=file[2]
    zmax=file[3]
    
    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0
    
    labeled, num_objects = ndimage.label(maxima)
    xy = np.array(ndimage.center_of_mass(data, labeled, range(1, num_objects+1)))
    
    plt.figure(100)
    plt.imshow(data)
    plt.savefig('./data.pdf', bbox_inches = 'tight')
    plt.autoscale(False)
    plt.title('Graph for maxima')
    plt.xlabel('Number of x points')
    plt.ylabel('Number of y points')
    plt.plot(xy[:, 1], xy[:, 0], 'ro')
    plt.savefig('./maxima.pdf', bbox_inches = 'tight')  
    plt.figure(count_del)
    plt.title('Final Plot')
    plt.xlabel('E (MeV)')
    plt.ylabel(r'$\Gamma$ (MeV)')
    plt.contourf(dataog[0],dataog[1], cut(dataog[2],graphcut/zmax))
    '''nothing bellow this is defined inside the function'''
    if masstf is True:    
        print("\n>>> channel \t thresholds")
        Mth=[i[0].real + i [1].real for i in MASS]
        for i in range(len(Mth)):
            if Mth[i] < Emax and Mth[i]> Emin:
                plt.axvline(x=Mth[i],color='r')
            print(i+1,'\t\t',round(Mth[i],3))  
    plt.savefig('./3Dgraph.pdf')
    return xy


#    
'''
Parameters
'''
from matrix import masses

graphcut=100
su3=False
cutoff=630
yfixed=1 # mass values
y=yfixed
xfixed=0. #type of renorm
x=xfixed
J=0.5
masstf=True
Emin=1300
Emax=2000
Ep=50 #number of points in energy
Wmin=0
Wmax=150
Wp=30 #number of points in width
ch1=None
ch2=ch1
RS=[[0 for i in range(dim(J))],[1 for i in range(dim(J))]]
RS2=[[0 for i in range(dim(J))],[1 for i in range(dim(J))]]
#RS[1][0]=0
'''
End of parameters
'''
XKcent=[]
MKcent=[]
WKcent=[]

count_del=0
while y<=yfixed:#
    XKcent.append(str(round(x,3)))
    MKcent.append([])
    WKcent.append([])    
    
    count_del+=1
    MASS=[masses(y,i,J) for i in range(4)]
    t1= t.time()
    useclusters=False
    nclust=6

    
    neighborhood_size = 5
    
    threshold = 0.00003
    
    
    dataog = pcalc.getpoles(Emin,Emax,Wmin,Wmax,Ep,Wp,ch1,ch2,x,cutoff,y,J,LnSRS=RS)
    if config.message != set():
        print(config.message)
    xy=plotmax(dataog,count_del,masstf,graphcut)

    
    if useclusters==True:
    #input("How many poles do you expect?\n")
        clusters = KMeans(n_clusters=nclust,random_state=0).fit(xy)
        Kcent=clusters.cluster_centers_
        plt.autoscale(False)
        plt.plot(Kcent[:, 1], Kcent[:, 0],'go')
        plt.savefig('./cluster_centers.png', bbox_inches = 'tight')
    #clusters.cluster_centers_
    else:
        Kcent=xy
        

    Kcent=[list(x) for x in list(Kcent)] #order list in terms of energy
    Kcent.sort(key=lambda L: L[1])
    print("\n>>>Fom graph:")
    for i in Kcent:
        print("Pole Position: "+str(i[1]*(Emax-Emin)/Ep+Emin)+"\t"+str(i[0]*(Wmax-Wmin)/Wp+Wmin))



#    import config#
    params=(x,cutoff,y,0.5)
    print("\n>>>Minimized:")
    POLES=[]
    for i in Kcent:
        p=pole(i[1]*(Emax-Emin)/Ep+Emin,i[0]*(Wmax-Wmin)/Wp+Wmin,params,vecRS=RS)
        with HiddenPrints():
            p.find()
        POLES.append(p)
        POLES.sort(key=lambda L: L.m)
        print("Pole Position: "+str(p.m)+"\t"+str(p.w))
    y+=1#0.493#

if len(POLES) > 3:    
    POLES=POLES[1:]
y-=1
print("\n>>>Couplings:")
if su3 is True:
    print("su(3): \t\t 1 \t\t 8 \t\t 8 \t\t 27")
    GABS=[3,2,0,1]
else:
    GABS=[0,1,2,3]
for i in POLES: 
    i.getcoupling(su3=su3)
    print(int(i.m),[i.gabs[x] for x in GABS])
    
    
t2= t.time()    
print("\nThe program took "+str(round(t2-t1,4))+" seconds to run.")
#    
#from matrix import mpi, mk, meta, fpi, fk, feta
#
#
##mp vs pole position
#print("\n",mpi[y],"\t",fpi[y],"\t",fk[y],"\t",feta[y],"\t",POLES[0].m,"\t",POLES[0].w,"\t",POLES[1].m,"\t", \
#  POLES[1].w,"\t",POLES[2].m,"\t", \
#  POLES[2].w)
#
##mpi vs su3 components 1 8 8_1
#if su3 is True:
#    print("\n",mpi[y],"\t",POLES[0].gabs[3],"\t",POLES[0].gabs[2],"\t",POLES[0].gabs[0],"\t",POLES[0].gabs[1], \
#      "\t",POLES[1].gabs[3],"\t", POLES[1].gabs[2],"\t",POLES[1].gabs[0],"\t",POLES[1].gabs[1],"\t", \
#      POLES[2].gabs[3],"\t",POLES[2].gabs[2],"\t",POLES[2].gabs[0],"\t",POLES[2].gabs[1])   

