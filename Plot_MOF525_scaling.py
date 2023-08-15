from __future__ import print_function
from ase import Atoms, Atom
from ase.optimize import BFGS
from ase.optimize import QuasiNewton
from ase.io import read
from ase.calculators.emt import EMT
from ase.visualize import view
from ase.io import write
import matplotlib.pyplot as plt
import numpy as np
from ase.build import molecule
import ase.db
from ase.db import connect
import os
from scipy import stats

###############################################################################
###   molecules
##############################################################################
db = ase.db.connect('data/molecules.db')
EH=db.get(mol='H2').energy*0.5
ECO=db.get(mol='CO').energy
EH2O=db.get(mol='H2O').energy


###############################################################################
###   MOF525
##############################################################################
folder='data/'
intermediates=['H','O','COOH','OH','OOH']
con = ase.db.connect(folder+'MOF525_3M_PW.db') # getting the slab energetics.

data_dict={}

size1=24;size2=20; size3=size2-10; sdot=150
cm = plt.cm.get_cmap('RdYlBu')

for Intermediates in intermediates:
        name_db=folder+'MOF525_3M_PW_%s.db' % (Intermediates)
        print(name_db)
        con1 = ase.db.connect(name_db)



        label_con1=[]
        Energy_con1=[]
        BG_con1=[]
        Metal_con1=[]
        

        
        for row in con1.select(relax='PW'):
            ESlab=con.get(metal=row.metal).energy
            #label_con1.append(row.ads)
            if Intermediates=='H':
                Energy_con1.append(row.energy-ESlab-EH)
            elif Intermediates=='O':
                Energy_con1.append(row.energy-ESlab+2*EH-EH2O)
            elif Intermediates=='OH':
                Energy_con1.append(row.energy-ESlab+EH-EH2O)
            elif Intermediates=='OOH':
                Energy_con1.append(row.energy-ESlab+3*EH-2*EH2O)
            elif Intermediates=='COOH':
                Energy_con1.append(row.energy-ESlab-ECO-EH2O+EH)
                
            BG_con1.append(row.bandgap)
            Metal_con1.append(row.metal)


        
        # plt.figure()
        # plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
        # plt.annotate(Intermediates, xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
        #              bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
        #              horizontalalignment='left', verticalalignment='bottom')
        # plt.scatter(Energy_con1,BG_con1, s=sdot,color='b')
        # for k in range(0,len(Energy_con1)):
        #     plt.text(Energy_con1[k],BG_con1[k], Metal_con1[k])

        # plt.xticks(fontsize=size2)
        # plt.yticks(fontsize=size2)
        # plt.xlabel('Energy [eV]', fontsize=size1)
        # plt.ylabel('Direct Bandgap [eV]', fontsize=size1)
        # #plt.title(Iron+'Fe with ' +Intermediates)
        # plt.ylim([0,2.0])
        # #plt.axis([-1, 1, -1, 10])
        # plt.savefig('MOF525_Metals_'+Intermediates+'.png', dpi=400, bbox_inches='tight')
        # plt.show()
        
# Save values to dictunary
        #data_dict[Intermediates+'_Label']=label_con1
        data_dict[Intermediates+'_Energy']=Energy_con1
        data_dict[Intermediates+'_Bandgap']=BG_con1
        data_dict[Intermediates+'_Metal']=Metal_con1

        
# #Check structure
# a = read(name_db+'@relax=PW,id=%s' % 109)
# val=-2
# chemsym=a[0].get_chemical_symbols()
# print(chemsym[val])
# chemsym[val]='S'
# a[0].set_chemical_symbols(chemsym)
# view(a[0])
        
###############################################################################
###   COOH vs H
##############################################################################
metals=['Co', 'Fe', 'Pd', 'Ni', 'Pt', 'Ir', 'Cu','Zn','Mn','Rh']

plt.figure()
plt.text(-1.1, 1.5, r'$\bf{b}$', ha='left', fontsize=size1)
x=[]
y=[]
for metal in metals:
        index_x=data_dict['H_Metal'].index(metal)
        index_y=data_dict['COOH_Metal'].index(metal)
        plt.scatter(data_dict['H_Energy'][index_x],data_dict['COOH_Energy'][index_y],color='b')
        plt.text(data_dict['H_Energy'][index_x],data_dict['COOH_Energy'][index_y],data_dict['COOH_Metal'][index_y],fontsize=size3)
        

        x.append(data_dict['H_Energy'][index_x])
        y.append(data_dict['COOH_Energy'][index_y])
        
gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
mn=np.min(x)
mx=np.max(x)
x1=np.linspace(mn,mx,500)
y1=gradient*x1+intercept
plt.plot(x1,y1,'-b')
plt.text(x1[-1]-1.5,y1[-1]-0.5,str(np.round(gradient,1))+'x'+str(np.round(intercept,1))+'\n $r^2$=' +str(np.round(r_value,1)),fontsize=size2-4)
            
plt.xticks(fontsize=size2)
plt.yticks(fontsize=size2)
plt.xlabel('$\Delta$E$_{H^*}$ [eV]',fontsize=size1)
plt.ylabel('$\Delta$E$_{COOH^*}$ [eV]',fontsize=size1)
plt.savefig('MOF525_Metals_COOH_vs_H.png', dpi=400, bbox_inches='tight')


###############################################################################
###   O vs OH
##############################################################################
plt.figure()
plt.text(0.25, 5.5, r'$\bf{c}$', ha='left', fontsize=size1)
x=[]
y=[]
z=[]
for metal in metals:
        index_x=data_dict['OH_Metal'].index(metal)
        index_y=data_dict['O_Metal'].index(metal)
        index_z=data_dict['OOH_Metal'].index(metal)
        plt.scatter(data_dict['OH_Energy'][index_x],data_dict['O_Energy'][index_y],color='b')
        plt.scatter(data_dict['OH_Energy'][index_x],data_dict['OOH_Energy'][index_z],color='k')
        plt.text(data_dict['OH_Energy'][index_x],data_dict['O_Energy'][index_y],data_dict['O_Metal'][index_y],fontsize=size3)
        plt.text(data_dict['OH_Energy'][index_x],data_dict['OOH_Energy'][index_z],data_dict['OOH_Metal'][index_z],fontsize=size3)        

        x.append(data_dict['OH_Energy'][index_x])
        y.append(data_dict['O_Energy'][index_y])
        z.append(data_dict['OOH_Energy'][index_z])
        
gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
mn=np.min(x)
mx=np.max(x)
x1=np.linspace(mn,mx,500)
y1=gradient*x1+intercept
plt.plot(x1,y1,'-b')
plt.text(x1[-1]-0.5,y1[-1]-2.0,str(np.round(gradient,1))+'x+'+str(np.round(intercept,1))+'\n $r^2$=' +str(np.round(r_value,1)), fontsize=size2-4, color='blue')

gradient, intercept, r_value, p_value, std_err = stats.linregress(x,z)
mn=np.min(x)
mx=np.max(x)
x1=np.linspace(mn,mx,500)
y1=gradient*x1+intercept
plt.plot(x1,y1,'-k')
plt.text(x1[-1]-1.5,y1[-1]-0.8,str(np.round(gradient,1))+'x+'+str(np.round(intercept,1))+'\n $r^2$=' +str(np.round(r_value,1)),fontsize=size2-4)

plt.xticks([0.5,1,1.5,2.0,2.5],[0.5,1.0,1.5,2.0,2.5],fontsize=size2)
plt.yticks([2.0,3,4,5],[2.0,3.0,4.0,5.0],fontsize=size2)
plt.xlabel('$\Delta$E$_{OH^*}$ [eV]',fontsize=size1)
plt.ylabel('{$\Delta$E$_{O^*}$, $\Delta$E$_{OOH^*}$}  [eV]',fontsize=size1-2)
plt.savefig('MOF525_Metals_O_vs_OH.png', dpi=400, bbox_inches='tight')


