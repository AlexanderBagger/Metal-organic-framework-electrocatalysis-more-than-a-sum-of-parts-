#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 09:50:09 2022

@author: kbv948
"""

from __future__ import print_function
from ase import Atoms
from ase.lattice import bulk
from ase.optimize.bfgs import BFGS
from ase.optimize import QuasiNewton
from ase.constraints import UnitCellFilter
from ase.constraints import FixAtoms, FixedPlane
import os
from ase.io import write, read
from ase.visualize import view
import numpy as np
from ase.io import Trajectory
from ase.lattice.surface import surface
from ase.dft.bee import BEEFEnsemble
import matplotlib.pyplot as plt
from numpy import genfromtxt

import os
import glob
####### Settings: ########

size1=24;size2=20; sdot=150



mylist=['NU-Por-4P-ftw_RPBE_bulk_PW500.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_1Fe.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_2Fe.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe.gpw.csv']
folder='data/DOS_data/'

mylist0=['NU-Por-4P-ftw_RPBE_bulk_PW500_1Fe_s2.gpw_SPIN0.csv',
        'NU-Por-4P-ftw_RPBE_bulk_PW500_2Fe_s01_s2.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN0.csv']
mylist1=['NU-Por-4P-ftw_RPBE_bulk_PW500_1Fe_s2.gpw_SPIN1.csv',
        'NU-Por-4P-ftw_RPBE_bulk_PW500_2Fe_s01_s2.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN1.csv']


names=['1Fe','2Fe','3Fe']


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#Density of States
plt.figure()
plt.figure(figsize=(7,7));

name1=folder+'NU-Por-4P-ftw_RPBE_bulk_PW500.gpw.csv'
data1=genfromtxt(name1,delimiter=',')
val11=find_nearest(data1[0,:],-1)
val12=find_nearest(data1[0,:],1)
data11=data1[:,val11:val12]
plt.plot(data11[0,:], (data11[1,:])/(data11[1,:]).max(),linewidth=4)
plt.text(-0.5,0.1,'MOF525', fontsize=size2)

for k in range(0,len(mylist0)):
        name1=folder+mylist0[k]
        data1=genfromtxt(name1,delimiter=',')
        name2=folder+mylist1[k]
        data2=genfromtxt(name2,delimiter=',')
        
        val11=find_nearest(data1[0,:],-1)
        val12=find_nearest(data1[0,:],1)
        #val21=find_nearest(data2[0,:],-1)
        #val22=find_nearest(data2[0,:],1)
        
        data11=data1[:,val11:val12]
        data12=data2[:,val11:val12]
        
        plt.plot(data11[0,:], (data11[1,:]+data12[1,:])/(data11[1,:]+data12[1,:]).max()+1*(k+1),linewidth=4)
        dx=0.0
        if names[k]=='1Fe':
            dx=-0.15
        elif names[k]=='2Fe' or names[k]=='3Fe':
            dx=-0.25            
        plt.text(-0.5+dx,1*(k+1)+0.1,names[k], fontsize=size2)

size1=28;size2=20; sdot=150
plt.xticks([-1,-0.5,0,0.5,1.0],['-1','-0.5','0.0','0.5','1.0'],fontsize=size2)
plt.yticks(fontsize=size2)
plt.plot([0,0],[-10,20],'k--')
plt.ylabel('DOS', fontsize=size1)
plt.xlabel('Energy [eV]', fontsize=size1)
#plt.legend(loc=1,fontsize=size2)
plt.axis([-1, 1, -0.1, 4.1])
plt.savefig('Plot_DOS_WithSpin.png', dpi=600, bbox_inches='tight')


mylist0=['NU-Por-4P-ftw_RPBE_bulk_PW500_3Cr_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Mn_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN0.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Co_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Ni_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Cu_s1_s1_s1.gpw_SPIN0.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Zn_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Rh_s1_s1_s1.gpw_SPIN0.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Pd_s1_s1_s1.gpw_SPIN0.csv', 'NU-Por-4P-ftw_RPBE_bulk_PW500_3Ir_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Pt_s1_s1_s1.gpw_SPIN0.csv']


mylist1=['NU-Por-4P-ftw_RPBE_bulk_PW500_3Cr_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Mn_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN1.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Co_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Ni_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Cu_s1_s1_s1.gpw_SPIN1.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Zn_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Rh_s1_s1_s1.gpw_SPIN1.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Pd_s1_s1_s1.gpw_SPIN1.csv', 'NU-Por-4P-ftw_RPBE_bulk_PW500_3Ir_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Pt_s1_s1_s1.gpw_SPIN1.csv']

mylist0=['NU-Por-4P-ftw_RPBE_bulk_PW500_3Mn_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN0.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Co_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Ni_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Cu_s1_s1_s1.gpw_SPIN0.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Zn_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Rh_s1_s1_s1.gpw_SPIN0.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Pd_s1_s1_s1.gpw_SPIN0.csv', 'NU-Por-4P-ftw_RPBE_bulk_PW500_3Ir_s1_s1_s1.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Pt_s1_s1_s1.gpw_SPIN0.csv']


mylist1=['NU-Por-4P-ftw_RPBE_bulk_PW500_3Mn_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN1.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Co_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Ni_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Cu_s1_s1_s1.gpw_SPIN1.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Zn_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Rh_s1_s1_s1.gpw_SPIN1.csv',
         'NU-Por-4P-ftw_RPBE_bulk_PW500_3Pd_s1_s1_s1.gpw_SPIN1.csv', 'NU-Por-4P-ftw_RPBE_bulk_PW500_3Ir_s1_s1_s1.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Pt_s1_s1_s1.gpw_SPIN1.csv']


names=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Density of States
plt.figure()
plt.figure(figsize=(7,8));
plt.text(-1.25, 10.2, r'$\bf{a}$', ha='left', fontsize=size1)
for k in range(0,len(mylist0)):
        name1=folder+mylist0[k]
        data1=genfromtxt(name1,delimiter=',')
        name2=folder+mylist1[k]
        data2=genfromtxt(name2,delimiter=',')
        
        val11=find_nearest(data1[0,:],-1)
        val12=find_nearest(data1[0,:],1)
        #val21=find_nearest(data2[0,:],-1)
        #val22=find_nearest(data2[0,:],1)
        
        data11=data1[:,val11:val12]
        data12=data2[:,val11:val12]
        
        plt.plot(data11[0,:], (data11[1,:]+data12[1,:])/(data11[1,:]+data12[1,:]).max()+1*k,linewidth=4)
        dx=0.0
        if names[k]=='Cr' or names[k]=='Zn':
            dx=0.15
        elif names[k]=='Mn' or names[k]=='Fe':
            dx=-0.2
        elif names[k]=='Co' or names[k]=='Ni' or names[k]=='Cu':
            dx=-0.1
        elif names[k]=='Ru' or names[k]=='Rh' or names[k]=='Pd' or names[k]=='Ir' or names[k]=='Pt':
            dx=0.9
            
        plt.text(-0.5+dx,1*k+0.1,names[k], fontsize=size2)



size1=28;size2=20; sdot=150
plt.xticks([-1,-0.5,0,0.5,1.0],['-1','-0.5','0.0','0.5','1.0'],fontsize=size2)
plt.yticks([],fontsize=size2)
plt.plot([0,0],[-10,20],'k--')
plt.ylabel('DOS', fontsize=size1)
plt.xlabel('Energy [eV]', fontsize=size1)
#plt.legend(loc=1,fontsize=size2)
plt.axis([-1, 1, -0.1, 10.1])
plt.savefig('Plot_DOS_Metals_Spin.png', dpi=600, bbox_inches='tight')


