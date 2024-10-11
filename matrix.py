# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 18:46:33 2018

@author: rpavao
"""

from cmath import phase,sqrt, exp,log, pi
import config



#==============================================================================
#Dm=[[-4.,sqrt(3./2.),-sqrt(1./2.)],[sqrt(3./2),-3.,-sqrt(27.0)],[-sqrt(1./2),-sqrt(27.0), -9.]]


'''
HERE J=1/2
'''

Dm12=[[-3,sqrt(3./2),-3./sqrt(2),0.],[sqrt(3./2),-4.,0.,-sqrt(3./2)], \
[-3./sqrt(2),0.,0.,3./sqrt(2)],[0.,-sqrt(3./2),3./sqrt(2),-3]]


Uma=[[0 for j in range(len(Dm12[i]))] for i in range(len(Dm12))]

Uma[1][1]=-0.5/sqrt(10.)
Uma[1][0]=-0.5*sqrt(3./5.)
Uma[1][2]=1.5*sqrt(0.3)
Uma[1][3]=0.5*sqrt(3./5.)
Uma[0][1]=-sqrt(3./5.)
Uma[0][0]=-1./sqrt(10.)
Uma[0][2]=-1./sqrt(5.)
Uma[0][3]=1./sqrt(10.)
Uma[2][1]=0.
Uma[2][0]=1./sqrt(2.)
Uma[2][2]=0.
Uma[2][3]=1./sqrt(2.)
Uma[3][1]=0.5*sqrt(1.5)
Uma[3][0]=-0.5
Uma[3][2]=-0.5/sqrt(2.)
Uma[3][3]=0.5

#
#Dm12=[[-3.,sqrt(1.5),-3./sqrt(2.),0.], \
#       [sqrt(1.5),-4.,0.,-sqrt(1.5)], \
#  [-3./sqrt(2.),0.,0., 3./sqrt(2.)], \
#  [0.,-sqrt(1.5), 3./sqrt(2.), -3.]]

# #MmB->Mass of Baryons, Mmm-> Mass of Mesons, Mf-> Meson decay constants

#units=4.63/10000.
def units(l):
    return (l/32.)/197

#our masses
L=[1,3.04,3.08,3.23,3.27]

'''japan mass '''
#config.message.add("japan mass")
#mpi=[x/units for x in [0.06391,0.1323,0.1895,0.2635,0.3220]]
#mk=[x/units for x in [0.2295,0.2747,0.2948,0.33,0.3622]]
#mnuc=[x/units for x in [0.4347,0.518,0.560,0.649,0.715]]
#mlamb=[x/units for x in [0.5165,0.593,0.643,0.695,0.762]]
#msig=[x/units for x in [0.5524,0.619,0.637,0.705,0.751]]

'''gubler mass with differnt L'''
config.message.add("gubler mass with differnt L")
mpi0=[138.039*units(1),0.1323,0.1895,0.2635,0.3220]
mk0=[495.3122*units(1),0.2747,0.2948,0.33,0.3622]
mnuc0=[938.919*units(1),0.518,0.560,0.649,0.715]
mlamb0=[1115.683*units(1),0.593,0.643,0.695,0.762]
msig0=[1193.1537*units(1),0.619,0.637,0.705,0.751]

mpi=[mpi0[i]/units(L[i]) for i in range(len(mpi0))]
mk=[mk0[i]/units(L[i]) for i in range(len(mpi0))]
mnuc=[mnuc0[i]/units(L[i]) for i in range(len(mpi0))]
mlamb=[mlamb0[i]/units(L[i]) for i in range(len(mpi0))]
msig=[msig0[i]/units(L[i]) for i in range(len(mpi0))]

meta = [sqrt((4.*mk[i]**2-mpi[i]**2)/3.) for i in range(len(mpi))]
meta[0]=547.862

mcas = [((msig[i]+3.*mlamb[i])*0.5-mnuc[i]) for i in range(len(mpi))]
mcas[0]=1318.285


'''f with doering parameters and L from doering and masses from philipp'''
config.message.add("f with doering parameters and L from doering and masses from philipp")
fpi=[92.3718,101.438,108.935,115.682,119.44]

fk=[112.592,117.275,119.654,121.388,122.362]

feta=[122.435,125.297,124.877,123.936,123.576]

'''f=90 '''
#config.message.add("f=90")
#fpi=[90 for i in range(5)]
#fk=[90 for i in range(5)]
#feta=[90 for i in range(5)]

'''D param with all init f0=90'''
#config.message.add("Relativistic Corrections")
##fpi=[89.9205,99.8836,107.763,115.294,118.763,120.86]
##fk=[103.275,113.874,119.486,122.04,122.09,120.86]
##feta=[102.027,105.269,111.919,118.369,120.881,120.86]

'''f with doering parameters'''
#config.message.add("f with doering parameters")
#fpi=[92.3704,101.997,109.832,117.092,120.296]
#fk=[112.61,118.012,120.545,122.58,122.982]
#feta=[122.409,126.336,125.964,125.216,124.182]





'''Doering masses and fpi'''
#config.message.add("Doering masses and fpi")
#mpi = [170.29, 282.84, 387.81, 515.56, 623.14]
#mk = [495.78, 523.26, 559.46, 609.75, 670.08]
#meta = [563.97, 581.72, 605.97, 638.07, 685.01]
#mnuc=[962.2,1058.7,1150.1,1274.5,1420.3]
#mlamb=[1135.8,1173.4,1261.0,1333.4,1434.2]
#msig=[1181.5,1235.5,1292.4,1353.5,1449.8]
#mcas=[1323.6,1332.8,1377.4,1401.8,1472.4]
#
#fpi=[94.5,102.5,109.5,116.3,120.1]
#fk=[113.2,116.1,118.5,120.6,121.9]
#feta=[122.1,122.3,122.6,122.4,122.6]

#
#
'''do not comment this'''
MmB=[mnuc,msig,mlamb,mcas]
Mmm=[mk,mpi,meta,mk]
Mf=[fk,fpi,feta,fk]

'''
HERE J=3/2
'''
Dm32=[]

#Dm32[5][10]=2./3.
#Dm32[10][5]=Dm32[5][10]

MmBc32=[]

Mmmc32=[]

Mfc32= []

MmBb32=[]

Mmmb32=[]

Mfb32= []


class Error_In_J(Exception):
    pass

def masses(ycb,i,J=0.5):
    if J==0.5:
        MmBi=MmB[i][ycb]
        Mmmi=Mmm[i][ycb]
        Mfi=Mf[i][ycb]
        return (MmBi,Mmmi,Mfi,Dm12)
        

#====================================================

def dim(J=0.5):
    if J==0.5:
        return len(Dm12)
    elif J==1.5:
        return len(Dm32)
    else:
        raise Error_In_J("You are using the wrong J")

#def nSRS(J):
#    SRS=[[0 for i in range(dim(J))],[1 for i in range(dim(J))]]
#    
##    RS=0
##    ch=1
##    if RS==0:
##        config.message.add('Changed ch'+str(ch+1)+' to SRS always.')
##    else:
##        config.message.add('Changed ch'+str(ch+1)+' to FRS always.')
##    SRS[RS][ch] = 1-RS
#    return SRS
 
#Col=map(sum, zip(*Dm))
#
#for i in range(len(Col)):
#    print i+1, Col[i].real
#    


def thresholds(ycb,J=0.5):
    for i in range(dim(J)):
        MASS=masses(ycb,i,J)
        print(i+1,"....",MASS[0]+MASS[1])
        
#def thresholds2(ycb,i,J):
#    TH_LIST=""
#    MASS=masses(ycb,i,J)
#    TH_LIST+=str(MASS[0]+MASS[1])+"\t"
#    return str(ycb)+"\t"+TH_LIST
#
#f=open("thresholds_mc_dependence",'w')
#
#for c in range(26):
#    y=0
#    while y<=1.:
#        f.write(thresholds2(y,c,1.5)+"\n")
#        y+=0.1
#    f.write("\n\n")
#
#f.close()
#    
