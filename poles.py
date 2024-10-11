#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 11:36:33 2019

@author: rpavao
"""
from cmath import phase,sqrt, exp,log, pi
import numpy as np
import numpy.linalg as nl
import scipy
import scipy.optimize as so
import scipy.integrate as sint
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from importlib import reload
import os, sys
import time as t
import config
from matrix import dim
from polescalc import T2,Tmat,getpoles
import matplotlib.pyplot as plt 
import sklearn 
from sklearn import cluster
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans
import importlib as ilib   
from plot_functions import get_xyz_files

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

class pole:
    '''  
    mass of pole
    width of pole
    *par = other parameters
    '''
    def __init__(self,mass,width,*par,vecRS=None):
        self.m=mass
        self.w=width
        try:
            while len(par)==1:
                    par=par[0]
        except TypeError:
            pass
        self.par=par
        self.J=self.par[3]
        self.dim=dim(self.J)
        if vecRS is None:
            self.RS=[[0 for i in range(dim(self.J))],[1 for i in range(dim(self.J))]]
        else:
            self.RS=vecRS
       
    def find(self,ch1=None,ch2=None):
        '''
        Based on the last values of mass and width, it finds the optimized on.
        No input needed, just put the Tmat(z,*par) function in the polescalc.py program.
        '''
        from polescalc import T2,Tmat,getpoles
        def Tinv(zl):            
            z=(zl[0]-0.5j*zl[1])**2
            return 1./T2(z,ch1,ch2,self.par,LnSRS=self.RS).real
        print('-----------------------------------')
        print('minit=',self.m,'\t winit=',self.w)
        sol=so.fmin(Tinv,np.array([float(self.m),float(self.w)]))
        self.m=sol[0]
        self.w=sol[1]
        print(self.m,self.w)
    
    def loop(self,ch,z=None,LnSRS=None):
        '''
        This calculates the loop function
        '''
        from functions import G
        if z is None:
            z=(self.m-.5j*self.w)**2
        if LnSRS is None:
            LnSRS=self.RS
        return G(z,ch,self.par,LnSRS=LnSRS)
    
    def kernel(self,ch1,ch2,z=None):
        '''
        This calculated the kernel
        '''
        from functions import Vf
        if z is None:
            z=(self.m-.5j*self.w)**2
        return Vf(z,ch1,ch2,self.par)
        
    def getcoupling(self,choch=None,diminit=0,dimvar=None,su3=False):
        '''
        Calculates de the couplings of the resonance (from Tii).
        No input needed, just put the Tmat(z,*par) function in the polescalc.py program.
        '''
        from matrix import Uma
        import numpy as np
        def Trxy(th,a,b):
            zr=(self.m-self.w*0.5j)
            z=(zr+r*exp(1.j*th))**2.
            if su3 is False:           
                Tm=Tmat(z,self.par,LnSRS=self.RS)
            else:
                A=np.array(Tmat(z,self.par,LnSRS=self.RS))
                B=np.array(Uma)
                Tm=np.matmul(B,np.matmul(A,np.transpose(B)))
                np.ndarray.tolist(Tm)      
            fun=r*exp(1.j*th)*Tm[a][b]/(2.*pi)
            return fun.real
        
        def Tixy(th,a,b):
            zr=(self.m-self.w*0.5j)
            z=(zr+r*exp(1.j*th))**2.
            if su3 is False:
                Tm=Tmat(z,self.par,LnSRS=self.RS)
            else:
                A=np.array(Tmat(z,self.par,LnSRS=self.RS))
                B=np.array(Uma)
                Tm=np.matmul(B,np.matmul(A,np.transpose(B)))
                np.ndarray.tolist(Tm)    
            fun=r*exp(1.j*th)*Tm[a][b]/(2.*pi)
            return fun.imag
        
        zr=(self.m-self.w*0.5j)**2.
        r=0.01
        COUP=[]
        COUPABS=[]
        if dimvar is None:
            print("Warning: relative signs of couplings not consistent")
            dimvar=self.dim
        if choch is not None:
            rhoch=sint.quad(Trxy,0,2.*pi,args=(choch-1,choch-1))[0]+sint.quad(Tixy, \
                           0,2.*pi,args=(choch-1,choch-1))[0]*1.j
        for i in range(diminit,dimvar):
            if choch is None:
                rho=sint.quad(Trxy,0,2.*pi,args=(i,i))[0]+sint.quad(Tixy,0,2.*pi,args=(i,i))[0]*1.j
                COUP.append(sqrt(rho))
                COUPABS.append(abs(sqrt(rho)).real)
            else:
                rho=(sint.quad(Trxy,0,2.*pi,args=(i,choch-1))[0]+sint.quad(Tixy,0,2.*pi, \
                     args=(i,choch-1))[0]*1.j)/sqrt(rhoch)
                COUP.append(rho)
                COUPABS.append(abs(rho).real)
        self.g=COUP
        self.gabs=COUPABS

    
    def distance(self,other,dimvar=None):#work on this "metric"
        '''
        Calculates the distance between your pole and the "other" poles. It does this by 
        checking the physical distance
        (mass and width) and also the diference between the couplings.
        This is still a work in progress.
        '''
        dcoup=0
        if dimvar is None:
            dimvar=self.dim
        for i in range(dimvar):
            try:
                dcoup+=(abs(self.g[i])-abs(other.g[i]))**2
            except AttributeError:
                self.getcoupling(dimvar)
                dcoup+=(abs(self.g[i])-abs(other.g[i]))**2
        dspace=((self.m-other.m)**2/max(self.m,other.m)+(self.w-other.w)**2/max(self.w,other.w))
        max1=max(map(abs,self.g))
        max2=max(map(abs,other.g))
        print("I have temporarily turned dspace to 0")
        return sqrt(0*dspace+dcoup/(max1+max2)).real
        
    def plot(self,Emin=None,Emax=None,Wmin=None,Wmax=None,Ep=40,Wp=40,par=None,just_plot=True, \
             npar=None,ch1=None,ch2=None):
        """
        If just_plot = True, this method plots a function (related to the determinant) of Tmat, within a given
        area around the poles. 
        If just_plot=False, it returns the matrix used to produce the plot (but no plot).
        """
        import matplotlib.pyplot as plt
        if Emin is None:
            Emin=self.m-20
        if Emax is None:
            Emax=self.m+20
        if Wmin is None:
            Wmin=self.w
        if Wmax is None:
            Wmax=self.w+20
        if par is None:
            par=self.par

        dataog = getpoles(Emin,Emax,Wmin,Wmax,Ep,Wp,ch1,ch2,par,LnSRS=self.RS)
        data=dataog[2]
        
        if just_plot==True:
            plt.figure()
            plt.contour(dataog[0],dataog[1], dataog[2])
            if npar is None:
                plt.savefig('./contour.pdf', bbox_inches = 'tight')
            else:
                plt.savefig('./Data/par'+str(npar+1)+'='+str(self.par[npar])+'.jpg', bbox_inches = 'tight')
        else:
            return  data
    
    def order(self,other):
        def kroneker(x,y):
            if x==y:
                return 1
            else:
                return 0
        Dself={}
        Dother={}
        for i in range(len(self.gabs)):
            Dself[self.gabs[i]]=i
            Dother[other.gabs[i]]=i
        Sself=sorted(self.gabs)
        Sother=sorted(other.gabs)
        Indexself=[]
        Indexother=[]
        for i in range(len(Sself)):
            Indexself.append(Dself[Sself[i]])
            Indexother.append(Dother[Sother[i]])
        return sum([kroneker(Indexself[i],Indexother[i]) for i in range(len(Indexself))])/self.dim

    def follow(self,parmax,dpar=None,n=0,name='datafollow',wantgraph=False,ch1=None,ch2=None):
        if dpar is None:
            dpar=0.1*parmax
        xpar=self.par[n]
        newpar=self.par
        massnew=self.m
        widthnew=self.w
        f=open(name,'w')
        while True:
            if dpar > 0 and xpar > parmax+dpar:
                break
            if dpar < 0 and xpar < parmax:
                break
            newpar=list(newpar)
            newpar[n]=xpar
            newpar=tuple(newpar)
            ptemp=pole(massnew,widthnew,newpar,vecRS=self.RS)
            with HiddenPrints():
                ptemp.find(ch1,ch2)
            massnew=ptemp.m
            widthnew=ptemp.w
            if ch1 is not None:
                ptemp.getcoupling(diminit=ch1-1,dimvar=ch1)
    #            print("We add to .dat file the ch=12 coupling")
                print(massnew,widthnew,xpar,ptemp.gabs[0])
                f.write(str(massnew)+"\t"+str(widthnew)+"\t"+str(xpar)+"\t"+str(ptemp.gabs[0])+"\n")
            else:
                print(massnew,widthnew,xpar)
                f.write(str(massnew)+"\t"+str(widthnew)+"\t"+str(xpar)+"\n")
            f.flush()
            xpar+=dpar
        f.close()
        if wantgraph is True:
            XYZ=get_xyz_files(name)
            fig,ax=plt.subplots(3)
            ax[0].plot(XYZ[0],XYZ[1])
            ax[1].plot(XYZ[2],XYZ[0])
            ax[2].plot(XYZ[2],XYZ[1])
            fig.subplots_adjust(hspace=0.6)
            fig.savefig('latest_follow_picture.pdf',bbox_inches='tight')
            
    def detectpoles(self,Emin=None,Emax=None,Wmin=None,Wmax=None,Ep=50,Wp=50,dE=100,dW=50 \
                    ,par=None,useclusters=False,nclust=None,neighborhood_size = 5,\
                    threshold = 0.003,npar=None,ch1=None,ch2=None):
        if Emin is None:
            Emin=self.m-dE
        if Emax is None:
            Emax=self.m+dE              
        if Wmin is None:
            if self.w <=dW:
                Wmin=0
            else:
                Wmin=self.w-dW
        if Wmax is None:
            Wmax=self.w+dW
        if par is None:
            par=self.par
            
        dataog = getpoles(Emin,Emax,Wmin,Wmax,Ep,Wp,ch1,ch2,par,LnSRS=self.RS)
        data=dataog[2]
        
        data_max = filters.maximum_filter(data, neighborhood_size)
        maxima = (data == data_max)
        data_min = filters.minimum_filter(data, neighborhood_size)
        diff = ((data_max - data_min) > threshold)
        maxima[diff == 0] = 0
        
        labeled, num_objects = ndimage.label(maxima)
        xy = np.array(ndimage.center_of_mass(data, labeled, range(1, num_objects+1)))
        
        
        plt.imshow(data)
#        plt.savefig('./Data/data'+str(self.par)+'.pdf', bbox_inches = 'tight')
        
        plt.autoscale(False)
        plt.plot(xy[:, 1], xy[:, 0], 'ro')
        if npar is None:
            plt.savefig('./maxima.pdf', bbox_inches = 'tight')
        else:
            plt.savefig('./Data/maxima_par'+str(npar+1)+'='+str(self.par[npar])+'.jpg', bbox_inches = 'tight')
        
        if useclusters==True:
        #input("How many poles do you expect?\n")
            clusters = KMeans(n_clusters=nclust,random_state=0).fit(xy)
            Kcent=clusters.cluster_centers_
            plt.autoscale(False)
            plt.plot(Kcent[:, 1], Kcent[:, 0],'go')
            if npar is None:
                plt.savefig('./cluster_centers_.pdf', bbox_inches = 'tight')
            else:
                plt.savefig('./Data/cluster_centers_par'+str(npar+1)+'='+str(self.par[npar])+'.jpg', bbox_inches = 'tight')
        #clusters.cluster_centers_
        else:
            Kcent=xy
        POLES=[]
        for i in Kcent:
            POLES.append([i[1]*(Emax-Emin)/Ep+Emin,i[0]*(Wmax-Wmin)/Wp+Wmin])
        return POLES        

    def evalallpoles(self,Emin=None,Emax=None,Wmin=None,Wmax=None, \
                    Ep=50,Wp=50,dE=100,dW=50,par=None,useclusters=False,nclust=None,\
                    neighborhood_size = 5,threshold = 0.003,POLES=None,ch1=None,ch2=None):
        t1=t.time()
        if POLES is None:
            POLES=self.detectpoles(Emin,Emax,Wmin,Wmax,Ep,Wp,dE,dW,par,useclusters, \
                                   nclust,neighborhood_size,threshold,ch1,ch2)
        DO=[]
        for i in POLES:
            ctemp=pole(i[0],i[1],self.par)
            with HiddenPrints():
                ctemp.find()
            ctemp.getcoupling()
            DO.append([self.distance(ctemp),self.order(ctemp)])
        t2=t.time()
        print("time:",t2-t1)
        return DO
    
#    def followregion(self,parmax,dpar=None,n=0,name='data_follow_region'):
#        if dpar is None:
#            dpar=0.1*parmax
#        xpar=self.par[n]
#        newpar=self.par
#        massnew=self.m
#        widthnew=self.w
#        f=open(name,'w')
#        while xpar<=parmax+dpar:
#            newpar=list(newpar)
#            newpar[n]=xpar
#            newpar=tuple(newpar)
#            ptemp=pole(massnew,widthnew,newpar)
#            with HiddenPrints():
#                ptemp.find()
#            f.write(str(ptemp.m)+"\t"+str(ptemp.w)+"\t"+str(xpar)+"\n")
#            print(ptemp.m,ptemp.w,xpar)
#            POLES=ptemp.detectpoles()
#            DO=ptemp.evalallpoles(POLES=POLES)
#            for i in range(len(POLES)):
#                if DO[i] == min(sorted(DO,key=lambda x: x[0])):
#                    massnew=POLES[i][0]
#                    widthnew=POLES[i][1]
#                    break
#            xpar+=dpar
#        f.close() 

#c7=pole(2804.83,20.74,0,1090,0.)
#c7.find()
#c7.follow(1.,0.1,2,"follow_c7-coup")
#
#
#import config
#print(config.message)


''' Pole c7 seems to end in M=6035. But that should be
the c5. I am now comparing the evolution of M vs ycb vs g_4 to
see which label is the correct. c5 should be the correct one
by group symmetry. Which pole is c7?
'''