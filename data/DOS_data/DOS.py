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

mylist = [f for f in glob.glob("*.csv")]
print(mylist)


mylist=['NU-Por-4P-ftw_RPBE_bulk_PW500.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_1Fe.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_2Fe.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe.gpw.csv']


# Density of States
plt.figure()
for struc, name, color in zip(mylist,['0Fe', '1Fe','2Fe','3Fe'],['k','b','r','c']):
        name1=struc
        data=genfromtxt(name1,delimiter=',')
        plt.plot(data[0,:], data[1,:],linewidth=4,color=color, label=name)

size1=28;size2=20; sdot=150
plt.xticks([-1,-0.5,0,0.5,1.0],['-1','-0.5','0.0','0.5','1.0'],fontsize=size2)
plt.yticks(fontsize=size2)
plt.ylabel('DOS', fontsize=size1)
plt.xlabel('Energy [eV]', fontsize=size1)
plt.legend(fontsize=size2)
plt.axis([-1, 1, -1, 20])
plt.savefig('Plot_DOS_NoSpin.png', dpi=600, bbox_inches='tight')


mylist0=['NU-Por-4P-ftw_RPBE_bulk_PW500.gpw.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_1Fe_s2.gpw_SPIN0.csv',
        'NU-Por-4P-ftw_RPBE_bulk_PW500_2Fe_s01_s2.gpw_SPIN0.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN0.csv']
mylist1=['NU-Por-4P-ftw_RPBE_bulk_PW500_1Fe_s2.gpw_SPIN1.csv',
        'NU-Por-4P-ftw_RPBE_bulk_PW500_2Fe_s01_s2.gpw_SPIN1.csv','NU-Por-4P-ftw_RPBE_bulk_PW500_3Fe_s1_s1_s1.gpw_SPIN1.csv']


# Density of States
plt.figure()
for struc, name, color in zip(mylist0,['0Fe', '1Fe','2Fe','3Fe'],['k','b','r','c']):
        name1=struc
        data=genfromtxt(name1,delimiter=',')
        plt.plot(data[0,:], data[1,:],linewidth=4,color=color, label=name)

for struc, name, color in zip(mylist1,['1Fe','2Fe','3Fe'],['b','r','c']):
        name1=struc
        data=genfromtxt(name1,delimiter=',')
        plt.plot(data[0,:], data[1,:],linewidth=4,color=color, linestyle=':')



size1=28;size2=20; sdot=150
plt.xticks([-1,-0.5,0,0.5,1.0],['-1','-0.5','0.0','0.5','1.0'],fontsize=size2)
plt.yticks(fontsize=size2)
plt.ylabel('DOS', fontsize=size1)
plt.xlabel('Energy [eV]', fontsize=size1)
plt.legend(loc=1,fontsize=size2)
plt.axis([-1, 1, -1, 20])
plt.savefig('Plot_DOS_WithSpin.png', dpi=600, bbox_inches='tight')





