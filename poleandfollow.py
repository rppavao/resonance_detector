#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 13:14:17 2019

@author: rpavao
"""
from poles import pole, HiddenPrints
from matrix import dim

'''
CHARM SECTOR
'''

#x_dcut=0
#y_cb=0
#cutoff=1400
#params12=(x_dcut,cutoff,y_cb,0.5)
#params32=(x_dcut,cutoff,y_cb,1.5)
#
#c1=pole(2699.4,12.6,params12)
#c2=pole(2702.8,177.8,params12)
#c3=pole(2733.,2.2,params12)
#c4=pole(2734.3,0.,params32)
#c5=pole(2772.9,83.7,params12)
#c6=pole(2775.4,0.6,params12)
#c7=pole(2804.8,20.7,params12)
#c8=pole(2819.7,32.4,params32)
#c9=pole(2845.2,44.,params32)
#
#CHARM=[c1,c2,c3,c4,c5,c6,c7,c8,c9]

'''
BOTTOM SECTOR
'''

x_dcut=0
cutoff=1400
y_cb=1
params=(x_dcut,cutoff,y_cb,0.5)




p1=pole(5873.83,0.,params)
p2=pole(5940.8,35.5,params)
p3=pole(5880.7,0,params)


'''
-------------------------------------------------------------
-------------------------------------------------------------
-------------------------------------------------------------
#'''
chosen=

#n_mass=chosen.m
#n_width=chosen.w


columns=""

columns=columns.split()

n_mass=float(columns[0])
n_width=float(columns[1])
n_y=1.

cutoff=float(columns[2])
n_x=1
J=chosen.J
RS=chosen.RS#[[0 for i in range(dim(J))],[1 for i in range(dim(J))]]
#RS[0][0]=1



b_n=pole(n_mass,n_width,n_x,cutoff,n_y,J,vecRS=RS)
#b_n.find()
#b_n.getcoupling()
b_n.follow(1394.2400000000052,-0.00001,1,"12-14_c4",3,3)



import config
print(config.message)