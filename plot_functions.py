# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 12:14:37 2019

@author: rpavao
"""

def get_xyz_files(file_name):
    f=open(file_name,'r')
    FILE=f.readlines()
    f.close()
    XYZ=[[] for i in range(len(FILE[0].split()))]
    for i in FILE:
        try: 
            xy=i.split()
            for j in range(len(XYZ)):
                if xy[0][0] !="#":
                    XYZ[j].append(float(xy[j]))
        except IndexError:
            print("IndexError")
            print(i)
            continue
    return XYZ
    
def get_xyz_function(f,a,b,npoints):
    X=[]
    Y=[]
    for i in range(npoints+1):
        x = a+(b-a)*i/npoints
        X.append(x)
        Y.append(f(x))
    return (X,Y)
        
def change_scale(my_list,factor):   
    my_list[:]=[x*factor for x in my_list]