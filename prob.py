  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 19:55:14 2019

@author: rpavao
"""
from cmath import sqrt
from poles import pole,HiddenPrints
from functions import G
from plot_functions import get_xyz_files,change_scale
if config.message != set():
    print(config.message)
def Div(f,x):
    h=x*10**(-5)
    return (f(x+h) - f(x-h))/(2.*h)


def prob(P):
#    with HiddenPrints():
    P.find()
    P.getcoupling(su3=False)
    print(P.g)
    zr=(P.m-.5j*P.w)**2
    return [-P.g[i]**2 *2*sqrt(zr)* Div(lambda x: P.loop(i+1,z=x),zr) for i in range(len(P.g))]
        
J=0.5
RS=[[0 for i in range(dim(J))],[1 for i in range(dim(J))]]
RS2=[[0 for i in range(dim(J))],[1 for i in range(dim(J))]]
RS2[1][0]=0



name="/with_relat_corrections/chiral_w_pmass_Ldoer"
XYZ=get_xyz_files("./results/"+name)
X=XYZ[0]
Msing=XYZ[4]
Wsing=XYZ[5]
LsingRS=[RS,RS,RS,RS,RS]
M1405=XYZ[6]
W1405=XYZ[7]
L1405RS=[RS2,RS2,RS,RS,RS]
M1670=XYZ[8]
W1670=XYZ[9]
L1670RS=[RS,RS,RS,RS,RS]

#f=open("./results/flavor"+name+"_sing","w")
#g=open("./results/flavor"+name+"_1405","w")
#h=open("./results/flavor"+name+"_1670","w")
#
#for y in range(len(X)):
#    par=(0,630,y,0.5)
#    Psing=prob(pole(Msing[y],Wsing[y],par,vecRS=LsingRS[y]))
#    P1405=prob(pole(M1405[y],W1405[y],par,vecRS=L1405RS[y]))
#    P1670=prob(pole(M1670[y],W1670[y],par,vecRS=L1670RS[y]))
#    linesing=str(X[y])+"\t"+str(Psing[0].real)+"\t"+str(Psing[1].real)+"\t"+str(Psing[2].real)+"\t"+str(Psing[3].real)+"\n"
#    line1405=str(X[y])+"\t"+str(P1405[0].real)+"\t"+str(P1405[1].real)+"\t"+str(P1405[2].real)+"\t"+str(P1405[3].real)+"\n"
#    line1670=str(X[y])+"\t"+str(P1670[0].real)+"\t"+str(P1670[1].real)+"\t"+str(P1670[2].real)+"\t"+str(P1670[3].real)+"\n"
#    f.writelines(linesing)
#    g.writelines(line1405)
#    h.writelines(line1670)
#    
#f.close()
#g.close()
#h.close()
y=0
par=(0,630,y,0.5)
Mp=Msing[y]
Wp=Wsing[y]
RSp=LsingRS[y]
p=pole(Mp,Wp,par,vecRS=RSp)    
#Psing=prob(sing)
#print(sum(Psing))
p.find()
p.getcoupling(choch=0)
zr=(p.m-0.5j*p.w)**2
#X=sum([-p.g[i]**2*2*sqrt(zr)*Div(lambda x: p.loop(i+1,z=x),zr) for i in range(4)])
#Z=sum([sum([-p.g[i]*p.g[j]*2*sqrt(zr)*Div(lambda x: p.kernel(i+1,j+1,z=x),zr)*p.loop(i+1,zr)*p.loop(j+1,zr) for i in range(4)]) for j in range(4)])
#print("X=",X.real)
#print("Z=",Z.real)
#print("X+Z=",Z.real+X.real)

print("\n---------------------------------------\n")
X=sum([-p.g[i]**2*Div(lambda x: p.loop(i+1,z=x**2),sqrt(zr)) for i in range(4)])
Z=sum([sum([-p.g[i]*p.g[j]*Div(lambda x: p.kernel(i+1,j+1,z=x**2),sqrt(zr))*p.loop(i+1,zr) \
            *p.loop(j+1,zr) for i in range(4)]) for j in range(4)])
print("X=",X.real)
print("Z=",Z.real)
print("X+Z=",Z.real+X.real)



''' test with one pole - doering y=1'''

#
#y=1
#par=(0,650,y,0.5)
#Mp=1783.5419717773484      
#Wp=38.95505375037527
#RSp=RS
#p=pole(Mp,Wp,par,vecRS=RSp)    
##Psing=prob(sing)
##print(sum(Psing))
#p.find()
#p.getcoupling(choch=0)
#zr=(p.m-0.5j*p.w)**2
##X=sum([-p.g[i]**2*2*sqrt(zr)*Div(lambda x: p.loop(i+1,z=x),zr) for i in range(4)])
##Z=sum([sum([-p.g[i]*p.g[j]*2*sqrt(zr)*Div(lambda x: p.kernel(i+1,j+1,z=x),zr)*p.loop(i+1,zr)*p.loop(j+1,zr) for i in range(4)]) for j in range(4)])
##print("X=",X.real)
##print("Z=",Z.real)
##print("X+Z=",Z.real+X.real)
#
#print("\n---------------------------------------\n")
#X=sum([-p.g[i]**2*Div(lambda x: p.loop(i+1,z=x**2),sqrt(zr)) for i in range(4)])
#Z=sum([sum([-p.g[i]*p.g[j]*Div(lambda x: p.kernel(i+1,j+1,z=x**2),sqrt(zr))*p.loop(i+1,zr) \
#            *p.loop(j+1,zr) for i in range(4)]) for j in range(4)])
#print("X=",X.real)
#print("Z=",Z.real)
#print("X+Z=",Z.real+X.real)
#
#for i in p.g:
#    print (abs(i))