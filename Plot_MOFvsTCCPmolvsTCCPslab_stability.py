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
import pandas as pd

# From Fend
#ESlab=-1908.364500 # Need PW500 number

###############################################################################
###   molecules
##############################################################################
db = ase.db.connect('data/molecules.db')
EH2=db.get(mol='H2').energy
EH=db.get(mol='H2').energy*0.5
ECO=db.get(mol='CO').energy
EH2O=db.get(mol='H2O').energy
#EN2=db.get(mol='N2').energy

###############################################################################
###   MOF525
##############################################################################
folder='data/'
intermediates=['H','O','COOH','OH','OOH']
con = ase.db.connect(folder+'MOF525_3M_PW.db') # getting the slab energetics.

MOF525_dict={}

size1=24;size2=20; sdot=150
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
            print(row.metal + ': ' + str(ESlab))
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
            print(row.metal + ' ' + Intermediates)
# Save values to dictunary
        #data_dict[Intermediates+'_Label']=label_con1
        MOF525_dict[Intermediates+'_Energy']=Energy_con1
        MOF525_dict[Intermediates+'_Bandgap']=BG_con1
        MOF525_dict[Intermediates+'_Metal']=Metal_con1


###############################################################################
###   Porphyrin_slab
##############################################################################
#folder='../data/Porphyrin_slab_DB/'
folder='data/'
intermediates=['H','O','COOH','OH','OOH']
con = ase.db.connect(folder+'Porphyrin_slab.db') # getting the slab energetics.
porphyrin_slab_dict={}
cm = plt.cm.get_cmap('RdYlBu')
for Intermediates in intermediates:
        name_db=folder+'Por_RPBE_slab_PW500_spin_ads_r1_1_%s.db' % (Intermediates)
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
            print(row.metal + ' ' + Intermediates)
        
# Save values to dictunary
        #data_dict[Intermediates+'_Label']=label_con1
        porphyrin_slab_dict[Intermediates+'_Energy']=Energy_con1
        porphyrin_slab_dict[Intermediates+'_Bandgap']=BG_con1
        porphyrin_slab_dict[Intermediates+'_Metal']=Metal_con1


###############################################################################
###   Porphyrin_mol
##############################################################################
folder='../data/Porphyrin_mol_DB/'
folder='data/'
intermediates=['H','O','COOH','OH','OOH']
con = ase.db.connect(folder+'Por_RPBE_mol_spin.db') # getting the slab energetics.
porphyrin_mol_dict={}
cm = plt.cm.get_cmap('RdYlBu')
for Intermediates in intermediates:
        name_db=folder+'Por_RPBE_mol_spin_%s.db' % (Intermediates)
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
        
# Save values to dictunary
        #data_dict[Intermediates+'_Label']=label_con1
        porphyrin_mol_dict[Intermediates+'_Energy']=Energy_con1
        porphyrin_mol_dict[Intermediates+'_Bandgap']=BG_con1
        porphyrin_mol_dict[Intermediates+'_Metal']=Metal_con1

###############################################################################
###   comparison
##############################################################################

metals=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']

Zr_E=-14.004314/2.0 # PW=500
ZrO2_E=-104.324657/4.0 #PW=500

size1=28;size2=20; sdot=150
E_mol_Zr=[]
E_slab_Zr=[]
E_mol_ZrO2=[]
E_slab_ZrO2=[]

for name in metals:
    
    m=name
    name_db=folder+'Por_RPBE_mol_spin.db'
    a = read(name_db+'@metal='+m)
    sym=a[0].symbols
    sym=a[0].get_chemical_symbols()
    sym_m=sym.count(m)
    sym_N=sym.count('N')
    sym_C=sym.count('C')
    sym_O=sym.count('O')
    sym_H=sym.count('H')
    
    
    con = ase.db.connect(name_db)
    Emol=con.get(metal=m).energy
    
    name_db=folder+'Porphyrin_slab.db'
    a = read(name_db+'@metal='+m)
    sym=a[0].symbols
    sym=a[0].get_chemical_symbols()
    slab_sym_m=sym.count(m)
    slab_sym_N=sym.count('N')
    slab_sym_C=sym.count('C')
    slab_sym_O=sym.count('O')
    slab_sym_H=sym.count('H')
    
    con = ase.db.connect(name_db)
    Eslab=con.get(metal=m).energy
    
    
    name_db=folder+'MOF525_3M_PW.db'
    a = read(name_db+'@metal='+m)
    sym=a[0].symbols
    sym=a[0].get_chemical_symbols()
    MOF_sym_m=sym.count(m)
    MOF_sym_N=sym.count('N')
    MOF_sym_C=sym.count('C')
    MOF_sym_O=sym.count('O')
    MOF_sym_H=sym.count('H')
    MOF_sym_Zr=sym.count('Zr')
    
    
    con = ase.db.connect(name_db)
    Emof=con.get(metal=m).energy
    
    Dmol_Zr=Emol-Emof/3.0+(MOF_sym_O/3-sym_O)*(EH2O-EH2)+2*Zr_E+(MOF_sym_H/3-sym_H)*EH
    Dslab_Zr=Eslab-Emof/3.0+(MOF_sym_O/3-slab_sym_O)*(EH2O-EH2)+2*Zr_E+(MOF_sym_H/3-slab_sym_H)*EH
    print(m)
    print('Mol stability (Zr): ' + str(Dmol_Zr))
    print('Slab stability (Zr): ' + str(Dslab_Zr))
    
    Dmol_ZrO2=Emol-Emof/3.0+2*ZrO2_E+(MOF_sym_O/3-sym_O-2*2)*(EH2O-EH2)+(MOF_sym_H/3-sym_H)*EH
    Dslab_ZrO2=Eslab-Emof/3.0+2*ZrO2_E+(MOF_sym_O/3-slab_sym_O-2*2)*(EH2O-EH2)+(MOF_sym_H/3-slab_sym_H)*EH
    print('Mol stability (ZrO2): ' + str(Dmol_ZrO2))
    print('Slab stability (ZrO2): ' + str(Dslab_ZrO2))
    E_mol_Zr.append(Dmol_Zr)
    E_slab_Zr.append(Dslab_Zr)
    E_mol_ZrO2.append(Dmol_ZrO2)
    E_slab_ZrO2.append(Dslab_ZrO2)
    


metals=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']
xtik=[1,2,3,4,5,6,7,8,9,10]
plt.figure()
plt.bar(np.asarray(xtik)-0.4, E_mol_Zr, width=0.4,color='r')
plt.bar(xtik, E_slab_Zr,width=0.4,color='b')
plt.plot([0.0,13.0],[0.0,0.0],'k--')
plt.ylim([-0.2,12])
plt.xlim([0.3,13])
plt.plot([10.0,10.4],[10.8,10.8],'k')
plt.text(10.4, 10.78, r'Slab', ha='left', fontsize=size2-4)
plt.plot([9.5,10.2],[5.2,5.2],'k')
plt.text(10.2, 5.18, r'Molecule', ha='left', fontsize=size2-4)
plt.plot([10.5,11.2],[0.0,2.12],'k')
plt.text(10.2, 2.12, r'MOF-525', ha='left', fontsize=size2-4)

plt.xticks(xtik, metals, fontsize=size2-2)
plt.yticks(fontsize=size2)
plt.xlabel('Metal loading', fontsize=size1)
plt.ylabel(r'$\Delta$ E [eV]', fontsize=size1)
plt.annotate(r'$\bf{b}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
plt.savefig('mol_slab_MOF_stability_Zr.png', dpi=400, bbox_inches='tight')

metals=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']
xtik=[1,2,3,4,5,6,7,8,9,10]
plt.figure()
plt.bar(np.asarray(xtik)-0.4, E_mol_ZrO2, width=0.4,color='r')
plt.bar(xtik, E_slab_ZrO2,width=0.4,color='b')
plt.plot([0.0,13.0],[0.0,0.0],'k--')
plt.ylim([-0.2,1.1])
plt.xlim([0.3,13])
plt.plot([10.0,10.4],[0.8,0.8],'k')
plt.text(10.4, 0.78, r'Slab', color='b',ha='left', fontsize=size2-4)
plt.plot([9.5,10.2],[0.2,0.2],'k')
plt.text(10.2, 0.18, r'Molecule', color='r',ha='left', fontsize=size2-4)
plt.plot([9.5,10.2],[0.0,0.-0.1],'k')
plt.text(10.2, -0.12, r'MOF-525', ha='left', fontsize=size2-4)

plt.xticks(xtik, metals, fontsize=size2-2)
plt.yticks(fontsize=size2)
plt.xlabel('Metal loading', fontsize=size1)
plt.ylabel(r'$\Delta$E [eV]', fontsize=size1)
plt.annotate(r'$\bf{a}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
plt.savefig('mol_slab_MOF_stability_ZrO2.png', dpi=400, bbox_inches='tight')

###############################################################################
###   Pourbaix stability
##############################################################################

E_MOF525=-1908.582664 # this includes 6H at N atoms

metals=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']


def seperate_string_number(string):
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i

        previous_character = i

        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups

folder='data/'
dict_bulk={}

bcc_slab=ase.db.connect(folder+'bcc_bulk_spin.db')
fcc_slab=ase.db.connect(folder+'fcc_bulk_spin.db')
hcp_slab=ase.db.connect(folder+'hcp_bulk_spin.db')


bcc_name=[]
bcc_energy=[]
fcc_name=[]
fcc_energy=[]
hcp_name=[]
hcp_energy=[]

for row in bcc_slab.select():
    name=row.formula
    print(name)
    data=seperate_string_number(name)
    Atoms1=read(folder+'bcc_bulk_spin.db@id=%s' %row.id)[0]
    sym=Atoms1.get_chemical_symbols()
    bcc_name.append(sym[0])
    if sym[0]=='Mn':
        bcc_energy.append(row.energy/2.0)
    else:
        bcc_energy.append(row.energy)
    
for row in fcc_slab.select():
    name=row.formula
    print(name)
    data=seperate_string_number(name)
    Atoms1=read(folder+'fcc_bulk_spin.db@id=%s' %row.id)[0]
    sym=Atoms1.get_chemical_symbols()
    fcc_name.append(sym[0])
    fcc_energy.append(row.energy)

for row in hcp_slab.select():
    name=row.formula
    print(name)
    data=seperate_string_number(name)
    Atoms1=read(folder+'hcp_bulk_spin.db@id=%s' %row.id)[0]
    sym=Atoms1.get_chemical_symbols()
    hcp_name.append(sym[0])
    hcp_energy.append(row.energy/2.0)


# Setting list of metal energetics.
metals=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']
m_name=[bcc_name[1],bcc_name[0],hcp_name[1],fcc_name[2],fcc_name[-2],hcp_name[0],fcc_name[0],fcc_name[3],fcc_name[1],fcc_name[-1]]
m_energy=[bcc_energy[1],bcc_energy[0],hcp_energy[1],fcc_energy[2],fcc_energy[-2],hcp_energy[0],fcc_energy[0],fcc_energy[3],fcc_energy[1],fcc_energy[-1]]


# from standard electro potential (data_page)
Mn_Vshe=-1.185  #(Mn^2+ + 2e^- -> Mn(s))
Fe_Vshe=-0.44   #(Fe^2+ + 2e^- -> Fe(s))
Co_Vshe=-0.28   #(Co^2+ + 2e^- -> Co(s))
Rh_Vshe=0.76    #(from somethere else)
Ir_Vshe=1.0     # from somewhere else)
Ni_Vshe=-0.25  #(Ni^2+ + 2e^- -> Ni(s))
Pd_Vshe=0.915  #(Pd^2+ + 2e^- -> Pd(s))
Pt_Vshe=1.188  #(Pt^2+ + 2e^- -> Pt(s))
Cu_Vshe=0.337  #(Cu^2+ + 2e^- -> Cu(s))
Zn_Vshe=-0.762 #(Zn^2+ + 2e^- -> Zn(s)

Vshe_energy=[Mn_Vshe,Fe_Vshe,Co_Vshe,Ni_Vshe,Cu_Vshe,Zn_Vshe,Rh_Vshe,Pd_Vshe,Ir_Vshe,Pt_Vshe]

# getting Emof stability.
EMOF=[]
for name in metals:
    name_db=folder+'MOF525_3M_PW.db'
    con = ase.db.connect(name_db)
    Emof=con.get(metal=name).energy
    EMOF.append(Emof-E_MOF525+6*EH)
    
    
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
metals=['Mn','Fe','Co','Ni','Cu','Zn','Rh','Pd','Ir','Pt']
xtik=[1,2,3,4,5,6,7,8,9,10]
ax=plt.figure()
plt.bar(np.asarray(xtik)-0.4, (np.asarray(EMOF)-3*np.asarray(m_energy))/3.0, width=0.4,color='r')
plt.bar(np.asarray(xtik)+0.0, (np.asarray(EMOF)-3*np.asarray(m_energy)-3*np.asarray(Vshe_energy))/3.0, width=0.4,color='b')
plt.plot([0.0,13.0],[0.0,0.0],'k--')
#plt.ylim([-0.2,1.1])
plt.xlim([0.3,13])

plt.plot([10.0,10.4],[-1.0,-1.0],'k')
plt.text(10.45, -1.05, r'vs m$^{2+}$', color='b',ha='left', fontsize=size2-4)
plt.plot([8.6,9.2],[0.5,0.5],'k')
plt.text(9.3, 0.48, r'vs m(s)', color='r',ha='left', fontsize=size2-4)

plt.xticks(xtik, metals, fontsize=size2-2)
plt.yticks(fontsize=size2)
plt.annotate(r'$\bf{b}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
plt.xlabel('Metal loading', fontsize=size1)
plt.ylabel([],fontsize=size1)
plt.ylabel(r'$\Delta$E [eV]', fontsize=size1)


#plt.annotate(r'$\bf{a}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
plt.savefig('metal_MOF_stability.png', dpi=400, bbox_inches='tight')



