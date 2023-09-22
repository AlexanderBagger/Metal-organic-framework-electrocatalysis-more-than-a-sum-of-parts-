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
intermediates=['H','O','COOH','OH','remove_H']


iron=['3']
data_dict={}
size1=28;size2=20; sdot=150
        
dis_max_test=0.8
cm = plt.cm.get_cmap('RdYlBu')
for Intermediates in intermediates:
    for Iron in iron:

        name_db=folder+'MOF525_%sFe_%s_PW_spin.db' % (Iron,Intermediates)
        print(name_db)
        con1 = ase.db.connect(name_db)
        
        if Iron=='0':
            ESlab=-1908.582664
            original=read('data/NU-Por-4P-ftw_RPBE_bulk_PW800.traj')
            POS_O=original.get_positions()

            
        elif Iron=='3':
            ESlab=-1913.48066
            original=read('data/NU-Por-4P-ftw_RPBE_bulk_PW800_3Fe.traj')
            POS_O=original.get_positions()

        label_con1=[]
        Energy_con1=[]
        BG_con1=[]
        dist_max=[]
        Fe_dist=np.empty(int(Iron)); Fe_dist.fill(100)
        Fe_Energy=np.empty(int(Iron)); Fe_Energy.fill(0)
        Fe_Bandgap=np.empty(int(Iron)); Fe_Bandgap.fill(0)
        Fe_Id=np.empty(int(Iron)); Fe_Bandgap.fill(0)
        
        for row in con1.select(relax='PW'):

# Getting the distances    
            a = read(name_db+'@relax=PW,ads=%s' % row.ads)
            newpos=a[0].get_positions()
            if (Intermediates=='H' or Intermediates=='O'):
                distances=np.sum((POS_O-newpos[:-1])**2,1)**(1/2)
                if distances.max()< dis_max_test:
                    if Intermediates=='H':
                        Energy_con1.append(row.energy-ESlab-EH)
                        dist_max.append(distances.max())
                        label_con1.append(row.ads)
                        BG_con1.append(row.bandgap)
                    elif Intermediates=='O':
                        Energy_con1.append(row.energy-ESlab+2*EH-EH2O)
                        dist_max.append(distances.max())
                        label_con1.append(row.ads)
                        BG_con1.append(row.bandgap)
                    
            elif (Intermediates=='CO' or Intermediates=='OH'):
                distances=np.sum((POS_O-newpos[:-2])**2,1)**(1/2)
                if distances.max()< dis_max_test:
                    if Intermediates=='OH':
                        Energy_con1.append(row.energy-ESlab+EH-EH2O)
                        dist_max.append(distances.max())
                        label_con1.append(row.ads)
                        BG_con1.append(row.bandgap)
                    elif Intermediates=='CO':
                        Energy_con1.append(row.energy-ESlab-ECO)
                        dist_max.append(distances.max())
                        label_con1.append(row.ads)
                        BG_con1.append(row.bandgap)
                
            elif Intermediates=='COOH':
                distances=np.sum((POS_O-newpos[:-4])**2,1)**(1/2)
                if distances.max()< dis_max_test:
                    Energy_con1.append(row.energy-ESlab-ECO-EH2O+EH)
                    dist_max.append(distances.max())
                    label_con1.append(row.ads)
                    BG_con1.append(row.bandgap)
                    
            elif Intermediates=='remove_H':
                POS_O=original.get_positions()
                val=row.ads
                POS_O=np.delete(POS_O, val, 0)
                distances=np.sum((POS_O-newpos)**2,1)**(1/2)
                if distances.max()< dis_max_test:
                    Energy_con1.append(row.energy-ESlab+EH)
                    dist_max.append(distances.max())
                    label_con1.append(row.ads)
                    BG_con1.append(row.bandgap)
                
            #print(distances.max())


        
# Save values to dictunary
        data_dict[Iron+'_'+Intermediates+'_Label']=label_con1
        data_dict[Iron+'_'+Intermediates+'_Energy']=Energy_con1
        data_dict[Iron+'_'+Intermediates+'_Bandgap']=BG_con1
        data_dict[Iron+'_'+Intermediates+'_Dist_max']=dist_max


###
# Adding Fe binding sites from seperate databases
#####
Fe_intermediates=['O','H','COOH','OH']

for Intermediates in Fe_intermediates:
        name_db=folder+'MOF525_3M_PW_%s.db' % (Intermediates)
        con = ase.db.connect(name_db) # getting the slab energetics.
        print(name_db)
        
        if Intermediates=='H':
            Fe_Energy=con.get(metal='Fe').energy-ESlab-EH
        elif Intermediates=='O':
            Fe_Energy=con.get(metal='Fe').energy-ESlab+2*EH-EH2O
        elif Intermediates=='OH':
            Fe_Energy=con.get(metal='Fe').energy-ESlab+EH-EH2O
        elif Intermediates=='COOH':
            Fe_Energy=con.get(metal='Fe').energy-ESlab-ECO-EH2O+EH

        Fe_Bandgap=con.get(metal='Fe').bandgap
        Fe_Id=con.get(metal='Fe').id
        
        # Getting the distances
        original=read('data/NU-Por-4P-ftw_RPBE_bulk_PW800_3Fe.traj')
        POS_O=original.get_positions()
        a = read(name_db+'@metal=Fe')
        newpos=a[0].get_positions()
        if (Intermediates=='H' or Intermediates=='O'):
            distances=np.sum((POS_O-newpos[:-1])**2,1)**(1/2)
        elif (Intermediates=='CO' or Intermediates=='OH'):
            distances=np.sum((POS_O-newpos[:-2])**2,1)**(1/2)
        elif (Intermediates=='COOH' or Intermediates=='OH'):
            distances=np.sum((POS_O-newpos[:-4])**2,1)**(1/2)

        Fe_dist=distances.max()        

        data_dict[Iron+'_'+Intermediates+'_Fe_dist']=Fe_dist
        data_dict[Iron+'_'+Intermediates+'_Fe_Energy']=Fe_Energy
        data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']=Fe_Bandgap
        data_dict[Iron+'_'+Intermediates+'_Fe_Id']=Fe_Id
        
#Check structure
#folder='data/'
#name_db=folder+'MOF525_3Fe_OH_PW_spin.db'
#a = read(name_db+'@relax=PW,id=%s' % 51)
#val=-2
#chemsym=a[0].get_chemical_symbols()
#print(chemsym[val])
#chemsym[val]='S'
#a[0].set_chemical_symbols(chemsym)
#view(a[0])

#Matching binding sites:
data_dict_match={}
for Intermediates in zip(['H','OH'],['COOH','O']):
#for Intermediates in zip(['H','H'],['COOH','CO']):
    for Iron in iron:
        print('Intermediates: ' ,Intermediates)
        print('Fe: ', Iron)

        datax=[]
        datay=[]
        dataGab1=[]
        dataGab2=[]
        id1=[]
        id2=[]
        dist=[]
        for k in range(0,len(data_dict[Iron+'_'+Intermediates[0]+'_Label'])):
            xmin=10
            ymin=10
            zminGab1=10
            zminGab2=10
            id1min=1
            id2min=1
            distmin=2
            for p in range(0,len(data_dict[Iron+'_'+Intermediates[1]+'_Label'])):
                    a1 = read(folder+'MOF525_'+Iron+'Fe_'+Intermediates[0]+'_PW_spin.db@relax=PW,ads=%s' % data_dict[Iron+'_'+Intermediates[0]+'_Label'][k])
                    newpos1=a1[0].get_positions()
                    a2 = read(folder+'MOF525_'+Iron+'Fe_'+Intermediates[1]+'_PW_spin.db@relax=PW,ads=%s' % data_dict[Iron+'_'+Intermediates[1]+'_Label'][p])
                    newpos2=a2[0].get_positions()
                    
                    if Intermediates[0]=='H':
                        val1=-1
                    elif Intermediates[0]=='OH':
                        val1=-2
                    if Intermediates[1]=='COOH':
                        val2=-4
                    elif Intermediates[1]=='O':
                        val2=-1
                    elif Intermediates[1]=='CO':
                        val2=-2
                    elif Intermediates[1]=='OH':
                        val2=-2
                    distances1=np.sum((newpos1[val1]-newpos2[val2])**2)**(1/2)
                    #print(distances1)
                    if distances1<distmin:
                            xmin=data_dict[Iron+'_'+Intermediates[0]+'_Energy'][k]
                            ymin=data_dict[Iron+'_'+Intermediates[1]+'_Energy'][p]
                            zminGab1=data_dict[Iron+'_'+Intermediates[0]+'_Bandgap'][k]
                            zminGab2=data_dict[Iron+'_'+Intermediates[1]+'_Bandgap'][p]
                            id1min=data_dict[Iron+'_'+Intermediates[0]+'_Label'][k]
                            id2min=data_dict[Iron+'_'+Intermediates[1]+'_Label'][p]
                            distmin=distances1
                        
                        # datax.append(data_dict['1_H_Energy'][k])
                        # datay.append(data_dict['1_COOH_Energy'][p])
                        # dataz.append(data_dict['1_H_Bandgap'][k])
                        # id1.append(data_dict['1_H_Label'][k])
                        # id2.append(data_dict['1_COOH_Label'][p])
            datax.append(xmin)
            datay.append(ymin)
            dataGab1.append(zminGab1)
            dataGab2.append(zminGab2)
            id1.append(id1min)
            id2.append(id2min)
            dist.append(distmin)      

        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datax']=datax
        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datay']=datay
        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_dataGab1']=dataGab1
        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_dataGab2']=dataGab2
        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_id1']=id1
        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_id2']=id2
        data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_dist']=dist
        
        
        
# Plot scaling :
for Intermediates in zip(['H','OH'],['COOH','O']):
#for Intermediates in zip(['H'],['COOH']):
    for Iron in iron:
        plt.figure()
        if Intermediates[0]=='H':
            plt.text(-0.4, 2.3, r'$\bf{c}$ Fe-MOF-525', ha='left', fontsize=size1)
        else:
            plt.text(-0.4, 3.65, r'$\bf{d}$ Fe-MOF-525', ha='left', fontsize=size1)
        
        #plt.annotate(Iron+'Fe', xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
        #              bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
        #              horizontalalignment='left', verticalalignment='bottom')
        plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
        
        plt.scatter(data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datax'],data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datay']
                    ,c=data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_dataGab1'], vmin=0, vmax=1.5, s=sdot, cmap=cm,edgecolors= "black")
        
        plt.scatter(data_dict[Iron+'_'+Intermediates[0]+'_Fe_Energy'],data_dict[Iron+'_'+Intermediates[1]+'_Fe_Energy']
                    ,c=data_dict[Iron+'_'+Intermediates[0]+'_Fe_bandgap'], vmin=0, vmax=1.5, s=sdot, cmap=cm,edgecolors= "black")
        
        cbar=plt.colorbar(ticks=[0.0, 0.5, 1.0,1.5])
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(size2)
        cbar.set_label('Bandgap [eV]', fontsize=size1)

        plt.annotate("Fe", xy=(data_dict[Iron+'_'+Intermediates[0]+'_Fe_Energy']+0.02,data_dict[Iron+'_'+Intermediates[1]+'_Fe_Energy']+0.02), xytext=(data_dict[Iron+'_'+Intermediates[0]+'_Fe_Energy']+0.15,data_dict[Iron+'_'+Intermediates[1]+'_Fe_Energy']+0.04),fontsize=size2-4,arrowprops=dict(arrowstyle="->"))
        plt.plot([-5,5],[-5,5],'k--')
        plt.xticks(fontsize=size2)
        plt.yticks(fontsize=size2)
        plt.xlabel('$\Delta$E$_{%s^*}$ [eV]' % Intermediates[0],fontsize=size1)
        plt.ylabel('$\Delta$E$_{%s^*}$ [eV]' % Intermediates[1],fontsize=size1)
        
        if Intermediates[0]=='OH':
            plt.ylim([1.0,3.5])
            plt.xlim([0.0,2.5])
        else:
            plt.ylim([-0.3,2.2])
            plt.xlim([0.0,2.5])

        plt.savefig(Iron+'Fe_with_'+Intermediates[0]+'vs'+Intermediates[1]+'_BG_SPIN.png', dpi=400, bbox_inches='tight')
        plt.show()



plt.figure(figsize=(18, 6))
plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5) 
xdat=[]
ydat=[]
zdat=[]
for k in range(0,len(data_dict['3_H_Label'])):
    xdat.append(-(data_dict['3_H_Energy'][k]+0.2)) # added H delta G correction
    ydat.append(data_dict['3_H_Bandgap'][k])
    zdat.append(data_dict['3_H_Dist_max'][k])

plt.scatter(-(data_dict['3_H_Fe_Energy']+0.2),data_dict['3_H_Fe_bandgap'],c=data_dict['3_H_Fe_dist'], vmin=0, vmax=1, s=sdot,marker='o', cmap=cm,edgecolors= "black")
p1=plt.scatter(xdat,ydat,c=zdat, vmin=0, vmax=1, s=sdot,marker='o', cmap=cm,edgecolors= "black", label=r'(H$^+$ +e$^-$) $\rightarrow$ H$^*$')

cbar=plt.colorbar(ticks=[0.0, 0.25, 0.5, 0.75, 1.0])
cbar.set_label('Distortion [Ang]', fontsize=size1)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(size2)

xdat2=[]
ydat2=[]
zdat2=[]
for k in range(0,len(data_dict['3_remove_H_Label'])):
    xdat2.append(data_dict['3_remove_H_Energy'][k]+0.2)
    ydat2.append(data_dict['3_remove_H_Bandgap'][k])
    zdat2.append(data_dict['3_remove_H_Dist_max'][k])
p2=plt.scatter(xdat2,ydat2,c=zdat2, vmin=0, vmax=1, s=sdot,marker='s', cmap=cm,edgecolors= "black", label=r'H$^*$ $\rightarrow$ (H$^+$ +e$^-$)')

xdat3=[]
ydat3=[]
zdat3=[]
for k in range(0,len(data_dict['3_OH_Label'])):
    xdat3.append(data_dict['3_OH_Energy'][k]+0.35-0.3) # added OH Delta G correction
    ydat3.append(data_dict['3_OH_Bandgap'][k])
    zdat3.append(data_dict['3_OH_Dist_max'][k])
p3=plt.scatter(xdat3,ydat3,c=zdat3, vmin=0, vmax=1, s=sdot,marker='^', cmap=cm,edgecolors= "black", label=r'H$_2$O $\rightarrow$ $^*$OH+(H$^+$ +e$^-$)')
plt.scatter(data_dict['3_OH_Fe_Energy']+0.35-0.3,data_dict['3_OH_Fe_bandgap'],c=data_dict['3_OH_Fe_dist'], vmin=0, vmax=1, s=sdot,marker='^', cmap=cm,edgecolors= "black")

xdat4=[]
ydat4=[]
zdat4=[]
for k in range(0,len(data_dict['3_O_Label'])):
    xdat4.append(data_dict['3_O_Energy'][k]+ 0.05)
    ydat4.append(data_dict['3_O_Bandgap'][k])
    zdat4.append(data_dict['3_O_Dist_max'][k])
p4=plt.scatter(xdat4,ydat4,c=zdat4, vmin=0, vmax=1, s=sdot,marker='p', cmap=cm,edgecolors= "black",label=r'H$_2$O $\rightarrow$ $^*$O+2(H$^+$ +e$^-$)')
plt.scatter(data_dict['3_O_Fe_Energy']+0.05,data_dict['3_O_Fe_bandgap'],c=data_dict['3_O_Fe_dist'], vmin=0, vmax=1, s=sdot,marker='p', cmap=cm,edgecolors= "black")


Iron='3'; Intermediates='H';
plt.annotate("Fe",weight="bold", xy=(-1*(data_dict[Iron+'_'+Intermediates+'_Fe_Energy']+0.2),data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']), xytext=(-1*(data_dict[Iron+'_'+Intermediates+'_Fe_Energy']+0.2)+0.15,data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']+0.05),fontsize=size2-4,arrowprops=dict(arrowstyle="->"))

Iron='3'; Intermediates='OH';
plt.annotate("Fe",weight="bold", xy=(data_dict[Iron+'_'+Intermediates+'_Fe_Energy']+0.35-0.3,data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']), xytext=(data_dict[Iron+'_'+Intermediates+'_Fe_Energy']-0.2+0.35-0.3,data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']+0.0),fontsize=size2-4,arrowprops=dict(arrowstyle="->"))


Iron='3'; Intermediates='O';
plt.annotate("Fe",weight="bold", xy=(data_dict[Iron+'_'+Intermediates+'_Fe_Energy']+0.05,data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']), xytext=(data_dict[Iron+'_'+Intermediates+'_Fe_Energy']+0.05,data_dict[Iron+'_'+Intermediates+'_Fe_bandgap']-0.1),fontsize=size2-4,arrowprops=dict(arrowstyle="->"))


def sigmoid(x,mi, mx): return mi + (mx-mi)*(lambda t: (1+200**(-t+2.0))**(-1) )( (x-mi)/(mx-mi) )
x = np.linspace(-4,7,1000)
plt.plot(x+0.35-0.3, sigmoid(x, 0.0, 0.4),'b-', lw=3, alpha=0.5)
plt.fill_between(x+0.35-0.3,sigmoid(x, 0.0, 0.4),color='b',alpha=0.2)

def sigmoid(x,mi, mx): return mi + (mx-mi)*(lambda t: (1+200**(t+1.0))**(-1) )( (x-mi)/(mx-mi) )
x = np.linspace(-4,7,1000)
plt.plot(x+0.2, sigmoid(x, 0.0, 0.4),'r-', lw=3, alpha=0.5)
plt.fill_between(x+0.2,sigmoid(x, 0.0, 0.4),color='r',alpha=0.2)

#plt.legend(fontsize=size2, loc=2, ncol=2)

legend1=plt.legend([p1], [r'(H$^+$ +e$^-$) $\rightarrow$ H$^*$'], loc=2, fontsize=size2)
legend2=plt.legend([p2,p3,p4], [r'H$^*$ $\rightarrow$ (H$^+$ +e$^-$)',r'H$_2$O $\rightarrow$ $^*$OH+(H$^+$ +e$^-$)',r'H$_2$O $\rightarrow$ $^*$O+2(H$^+$ +e$^-$)'], loc=1, fontsize=size2)
plt.gca().add_artist(legend1)
plt.text(-1.8, 1.05, r'$\bf{b}$ Fe-MOF-525', ha='left', fontsize=size2+14)

plt.text(-1.48, 0.02, r'Reduction', ha='left', fontsize=size2)
plt.text(1.58, 0.02, r'Oxidation', ha='left', fontsize=size2)


plt.xticks(fontsize=size2)
plt.yticks(fontsize=size2)
plt.xlabel('E / V$_{RHE}$ [eV]', fontsize=size1)
plt.ylabel('Bandgap [eV]', fontsize=size1)
plt.plot([0,0],[0,1],'k-',linewidth=4)
#plt.title(Iron+'Fe with ' +Intermediates)
plt.ylim([0,1.0])
plt.xlim([-1.5,2.0])
#plt.axis([-1, 1, -1, 10])
#plt.savefig(Iron+'Fe_with_'+Intermediates+'_SPIN.png', dpi=400, bbox_inches='tight')
plt.savefig('3Fe_with_SPIN_potential_map.png', dpi=400, bbox_inches='tight')
plt.show()

####################################
### Viewing lowest lying structures
####################################

# Iron='3'
# # Finding structure
# ymin_arg=np.argmin(data_dict_match[Iron+'_'+'H'+'_'+'COOH'+'_datay'])
# # Finding id's
# xid=data_dict_match[Iron+'_'+'H'+'_'+'COOH'+'_id1'][ymin_arg]
# yid=data_dict_match[Iron+'_'+'H'+'_'+'COOH'+'_id2'][ymin_arg]

# name_db='data/MOF525_3Fe_H_PW_spin.db'
# # #Check structure
# a = read(name_db+'@ads=%s' % xid)
# val=-1
# chemsym=a[0].get_chemical_symbols()
# print(chemsym[val])
# chemsym[val]='S'
# a[0].set_chemical_symbols(chemsym)
# view(a[0])

# name_db='data/MOF525_3Fe_COOH_PW_spin.db'
# # #Check structure
# a = read(name_db+'@ads=%s' % yid)
# val=-4
# chemsym=a[0].get_chemical_symbols()
# print(chemsym[val])
# chemsym[val]='S'
# a[0].set_chemical_symbols(chemsym)
# view(a[0])

Iron='3'
# Finding structure
ymin_arg=np.argmin(data_dict_match[Iron+'_'+'OH'+'_'+'O'+'_datax'])
# Finding id's
xid=data_dict_match[Iron+'_'+'OH'+'_'+'O'+'_id1'][ymin_arg]
yid=data_dict_match[Iron+'_'+'OH'+'_'+'O'+'_id2'][ymin_arg]

name_db='data/MOF525_3Fe_OH_PW_spin.db'
# #Check structure
a = read(name_db+'@ads=%s' % xid)
val=-2
chemsym=a[0].get_chemical_symbols()
print(chemsym[val])
chemsym[val]='S'
a[0].set_chemical_symbols(chemsym)
view(a[0])

name_db='data/MOF525_3Fe_O_PW_spin.db'
# #Check structure
a = read(name_db+'@ads=%s' % yid)
val=-1
chemsym=a[0].get_chemical_symbols()
print(chemsym[val])
chemsym[val]='S'
a[0].set_chemical_symbols(chemsym)
view(a[0])


# Check special structure

ymin_arg=np.argmin(data_dict['3_O_Energy'])
id=data_dict['3_O_Label'][ymin_arg]

name_db='data/MOF525_3Fe_O_PW_spin.db'
# #Check structure
a = read(name_db+'@ads=%s' % id)
val=-1
chemsym=a[0].get_chemical_symbols()
print(chemsym[val])
chemsym[val]='S'
a[0].set_chemical_symbols(chemsym)
view(a[0])





