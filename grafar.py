#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 16:38:43 2019

@author: rpavao
"""
#%matplotlib qt
import matplotlib.pyplot as plt
from plot_functions import get_xyz_files,change_scale
from poles import pole, HiddenPrints
from matrix import dim
from math import sqrt

units=4.63/10000.

#XYZ=get_xyz_files("./results/chiral_w_pmass_Ldoer")
#X=XYZ[0]
#Lsing=XYZ[4]
#L1405=XYZ[6]
#L1670=XYZ[8]
#
#
#XYZ=get_xyz_files("./results/paper1405_p")
#Xgub=[sqrt(x)/units for x in XYZ[0]]
#G1405=XYZ[1]
#d1405=XYZ[2]
#
#XYZ=get_xyz_files("./results/paper1670_p")
#G1670=XYZ[1]
#d1670=XYZ[2]
#
#XYZ=get_xyz_files("./results/paper1800_p")
#G1800=XYZ[1]
#d1800=XYZ[2]
#
#
#change_scale(G1405,1/units)
#change_scale(d1405,1/units)
#change_scale(G1670,1/units)
#change_scale(d1670,1/units)
#change_scale(G1800,1/units)
#change_scale(d1800,1/units)
#
#XYZ=get_xyz_files("./results/adelaide_data")
#Xad=XYZ[0]
#A1405=XYZ[1]
#dA1405=XYZ[2]
#Xad2=XYZ[0][-2:]
#A1670=XYZ[3][-2:]
#dA1670=XYZ[4][-2:]
#
#
##plt.plot(X,Lsing,c='k')
###plt.scatter(X,Lsing)
##
##plt.plot(X,L1405,c='k')
###plt.scatter(X,L1405)
##
##plt.plot(X,L1670,c='k')
###plt.scatter(X,L1670)
##
##
##plt.errorbar(Xgub,G1405,yerr=d1405,fmt='o')
##
##plt.errorbar(Xgub,G1670,d1670,fmt='x')
##
##plt.errorbar(Xgub,G1800,d1800,fmt='s')
##
##
##plt.title(r"Japan data ($f_0=79.2$)")
##
###plt.legend(loc=1,ncol=1, fancybox=False, shadow=False, \
###           prop=fontP,frameon=False,handlelength=2)
##
##plt.xlabel(r'$m_{\pi}$ (MeV)')
##plt.ylabel(r'$M_{\Lambda}$ (MeV)')
###plt.savefig("./results/japan_data_doering.pdf",bbox_inches='tight')
##
#
#
#
##plt.plot(X,Lsing,c='k')
###plt.scatter(X,Lsing)
##
##plt.plot(X,L1405,c='k')
###plt.scatter(X,L1405)
##
##plt.plot(X,L1670,c='k')
###plt.scatter(X,L1670)
##
##
##plt.errorbar(Xad,A1405,dA1405,fmt='o',label=r"$\Lambda(1405)_{adelaide}$")
##
##plt.errorbar(Xad2,A1670,dA1670,fmt='x',label=r"$\Lambda(1670)_{adelaide}$")
#
##
##
##plt.title(r"Adelaide data")
##
###plt.legend(loc=1,ncol=1, fancybox=False, shadow=False, \
###           prop=fontP,frameon=False,handlelength=2)
##
##plt.xlabel(r'$m_{\pi}$ (MeV)')
##plt.ylabel(r'$M_{\Lambda}$ (MeV)')
#
##plt.savefig("./results/adelaide_data_chiral_massesdoering.pdf",bbox_inches='tight')
#
#from matplotlib.font_manager import FontProperties
#fontP = FontProperties()
#fontP.set_size(14)
#
#XYZ=get_xyz_files("./results/paper1405_p_Ldoer")
#Xgub=XYZ[0]
#G1405=XYZ[1]
#d1405=XYZ[2]
#
#XYZ=get_xyz_files("./results/paper1670_p_Ldoer")
#G1670=XYZ[1]
#d1670=XYZ[2]
#
#XYZ=get_xyz_files("./results/paper1800_p_Ldoer")
#G1800=XYZ[1]
#d1800=XYZ[2]
#
#
#plt.plot(X,Lsing,c='k')
##plt.scatter(X,Lsing)
#
#plt.plot(X,L1405,c='k')
##plt.scatter(X,L1405)
#
#plt.plot(X,L1670,c='k')
##plt.scatter(X,L1670)
#
#
#plt.errorbar(Xgub,G1405,yerr=d1405,fmt='o',label=r"$\Lambda(1405)_{japan}$")
#
#plt.errorbar(Xgub,G1670,d1670,fmt='x',label=r"$\Lambda(1670)_{japan}$")
##
#plt.errorbar(Xgub,G1800,d1800,fmt='s',label=r"$\Lambda(1800)_{japan}$")
#
##
##plt.errorbar(Xad,A1405,dA1405,fmt='o',label=r"$\Lambda(1405)_{adelaide}$")
##
##plt.errorbar(Xad2,A1670,dA1670,fmt='x',label=r"$\Lambda(1670)_{adelaide}$")
#
#
#plt.title(r"Japan data with L from Doering")
##
#plt.legend(loc=1,ncol=1, fancybox=False, shadow=False, \
#           prop=fontP,frameon=False,handlelength=2,bbox_to_anchor=(1.5, 1.))
##
#plt.xlabel(r'$m_{\pi}$ (MeV)')
#plt.ylabel(r'$M_{\Lambda}$ (MeV)')
#plt.savefig("./results/chiral_japan_fandLdoering.pdf",bbox_inches='tight')


''' flavor content '''

def norm(L,j):
    s=0
    for i in L:
       s+=i[j]
    return 1

XYZ=get_xyz_files("./results/flavor_chiral_gublermasses")
#XYZ=get_xyz_files("./results/flavor_chiral_Ldoering")
#XYZ=get_xyz_files("./results/flavor_everything_adelaide")
X=XYZ[0]

g1_sing=XYZ[1]
g8_sing=XYZ[2]
g8p_sing=XYZ[3]
g27_sing=XYZ[4]

g1_1405=XYZ[5]
g8_1405=XYZ[6]
g8p_1405=XYZ[7]
g27_1405=XYZ[8]

g1_1670=XYZ[9]
g8_1670=XYZ[10]
g8p_1670=XYZ[11]
g27_1670=XYZ[12]

P=3

LofL=[[g1_sing,g8_sing,g8p_sing,g27_sing], [g1_1405,g8_1405,g8p_1405,g27_1405], \
[g1_1670,g8_1670,g8p_1670,g27_1670]]


G1=[LofL[P-1][0][i]/norm(LofL[P-1],i) for i in range(len(LofL[P-1][0]))]
G8=[LofL[P-1][1][i]/norm(LofL[P-1],i) for i in range(len(LofL[P-1][1]))]
G8p=[LofL[P-1][2][i]/norm(LofL[P-1],i) for i in range(len(LofL[P-1][2]))]

plt.scatter(X,G1,label='singlet')

plt.scatter(X,G8,label='octet 1')

plt.scatter(X,G8p,label='octet 2')


plt.title("IRREPS")
#
plt.legend(loc=1,ncol=1, fancybox=False, shadow=False, \
           prop=fontP,frameon=False,handlelength=2,bbox_to_anchor=(1.5, 1.))
#
plt.xlabel(r'$m_{\pi}$ (MeV)')
plt.ylabel(r'$g_{irrep}$')
plt.savefig("./results/irreps.pdf",bbox_inches='tight')