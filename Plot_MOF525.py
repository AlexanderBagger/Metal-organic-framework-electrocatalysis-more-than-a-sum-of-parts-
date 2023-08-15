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
iron=['0']

data_dict={}

size1=28;size2=20; sdot=150
        
cm = plt.cm.get_cmap('RdYlBu')
for Intermediates in intermediates:
    for Iron in iron:

        name_db=folder+'MOF525_%sFe_%s_PW.db' % (Iron,Intermediates)
        print(name_db)
        con1 = ase.db.connect(name_db)
        
        if Iron=='0':
            ESlab=-1908.582664
            original=read('data/NU-Por-4P-ftw_RPBE_bulk_PW800.traj')
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
            label_con1.append(row.ads)
            
            if Intermediates=='H':
                Energy_con1.append(row.energy-ESlab-EH)
            elif Intermediates=='remove_H':
                Energy_con1.append(row.energy-ESlab+EH)
            elif Intermediates=='O':
                Energy_con1.append(row.energy-ESlab+2*EH-EH2O)
            elif Intermediates=='OH':
                Energy_con1.append(row.energy-ESlab+EH-EH2O)
            elif Intermediates=='CO':
                Energy_con1.append(row.energy-ESlab-ECO)
            elif Intermediates=='COOH':
                Energy_con1.append(row.energy-ESlab-ECO-EH2O+EH)
                
            BG_con1.append(row.bandgap)

# Getting the distances    
            a = read(name_db+'@relax=PW,ads=%s' % row.ads)
            newpos=a[0].get_positions()
            if (Intermediates=='H' or Intermediates=='O'):
                distances=np.sum((POS_O-newpos[:-1])**2,1)**(1/2)
            elif (Intermediates=='CO' or Intermediates=='OH'):
                distances=np.sum((POS_O-newpos[:-2])**2,1)**(1/2)
            elif (Intermediates=='COOH' or Intermediates=='OH'):
                distances=np.sum((POS_O-newpos[:-4])**2,1)**(1/2)
            elif Intermediates=='remove_H':
                POS_O=original.get_positions()
                val=row.ads
                POS_O=np.delete(POS_O, val, 0)
                distances=np.sum((POS_O-newpos)**2,1)**(1/2)
            #print(distances.max())
            dist_max.append(distances.max())


        #atomsobj=con1.get_atoms(ads=395,relax='PW')
        #chemsym=atomsobj.get_chemical_symbols()
        #chemsym[-1]='S'
        #atomsobj.set_chemical_symbols(chemsym)
        #view(atomsobj)


        
        # plt.figure()
        # plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
        # plt.annotate(Iron+'Fe with '+Intermediates, xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
        #              bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
        #              horizontalalignment='left', verticalalignment='bottom')
        # plt.scatter(Energy_con1,BG_con1,c=dist_max, vmin=0, vmax=1, s=sdot, cmap=cm,edgecolors= "black")

        # cbar=plt.colorbar()
        # cbar.set_label('Max dist [Ang]', fontsize=size1)
        # plt.xticks(fontsize=size2)
        # plt.yticks(fontsize=size2)
        # plt.xlabel('Energy [eV]', fontsize=size1)
        # plt.ylabel('Direct Bandgap [eV]', fontsize=size1)
        # #plt.title(Iron+'Fe with ' +Intermediates)
        # plt.ylim([0,1.5])
        # #plt.axis([-1, 1, -1, 10])
        # plt.savefig(Iron+'Fe_with_'+Intermediates+'.png', dpi=400, bbox_inches='tight')
        # plt.show()
        
# Save values to dictunary
        data_dict[Iron+'_'+Intermediates+'_Label']=label_con1
        data_dict[Iron+'_'+Intermediates+'_Energy']=Energy_con1
        data_dict[Iron+'_'+Intermediates+'_Bandgap']=BG_con1
        data_dict[Iron+'_'+Intermediates+'_Dist_max']=dist_max

        
# #Check structure
# a = read(name_db+'@relax=PW,id=%s' % 109)
# val=-2
# chemsym=a[0].get_chemical_symbols()
# print(chemsym[val])
# chemsym[val]='S'
# a[0].set_chemical_symbols(chemsym)
# view(a[0])

# # Matching binding sites:
data_dict_match={}
for Intermediates in zip(['H','OH'],['COOH','O']):
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
            distmin=10
            for p in range(0,len(data_dict[Iron+'_'+Intermediates[1]+'_Label'])):
                    a1 = read(folder+'MOF525_'+Iron+'Fe_'+Intermediates[0]+'_PW.db@relax=PW,ads=%s' % data_dict[Iron+'_'+Intermediates[0]+'_Label'][k])
                    newpos1=a1[0].get_positions()
                    a2 = read(folder+'MOF525_'+Iron+'Fe_'+Intermediates[1]+'_PW.db@relax=PW,ads=%s' % data_dict[Iron+'_'+Intermediates[1]+'_Label'][p])
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
    for Iron in iron:
        plt.figure()
        if Intermediates[0]=='H':
            plt.text(-1.0, 1.35, r'$\bf{a}$ MOF-525', ha='left', fontsize=size1)
        else:
            plt.text(1.0, 6.5, r'$\bf{b}$ MOF-525', ha='left', fontsize=size1)

        plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
        plt.scatter(data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datax'],data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datay']
                    ,c=data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_dataGab1'], vmin=0, vmax=1.5, s=sdot, cmap=cm,edgecolors= "black")
        cbar=plt.colorbar(ticks=[0.0, 0.5, 1.0,1.5])
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(size2)
        cbar.set_label('Bandgap [eV]', fontsize=size1-4)

        
        plt.xticks(fontsize=size2)
        if Intermediates[0]=='H':
            plt.yticks([-1,0,1],['-1','0','1'],fontsize=size2)
        else:
            plt.yticks(fontsize=size2)
        plt.xlabel('$\Delta$E$_{%s^*}$ [eV]' % Intermediates[0],fontsize=size1-4)
        plt.ylabel('$\Delta$E$_{%s^*}$ [eV]' % Intermediates[1],fontsize=size1-4)
        

        plt.savefig(Iron+'Fe_with_'+Intermediates[0]+'vs'+Intermediates[1]+'_BG.png', dpi=400, bbox_inches='tight')
        plt.show()
        
        plt.figure()
        plt.annotate(Iron+'Fe', xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                      bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                      horizontalalignment='left', verticalalignment='bottom')
        plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
        plt.scatter(data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datax'],data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_datay']
                    ,c=data_dict_match[Iron+'_'+Intermediates[0]+'_'+Intermediates[1]+'_dist'], vmin=0, vmax=1.5, s=sdot, cmap=cm,edgecolors= "black")
        cbar=plt.colorbar(ticks=[0.0, 0.5, 1.0,1.5])
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(size2)
        cbar.set_label('dist [Ang]', fontsize=size1)
        
        plt.xticks(fontsize=size2)
        if Intermediates[0]=='H':
            plt.yticks([-1,0,1],['-1','0','1'],fontsize=size2)
        else:
            plt.yticks(fontsize=size2)
        plt.xlabel('$\Delta$E$_{%s^*}$ [eV]' % Intermediates[0],fontsize=size1)
        plt.ylabel('$\Delta$E$_{%s^*}$ [eV]' % Intermediates[1],fontsize=size1)
        plt.savefig(Iron+'Fe_with_'+Intermediates[0]+'vs'+Intermediates[1]+'_dist.png', dpi=400, bbox_inches='tight')
        plt.show()




plt.figure(figsize=(18, 6))
plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5) 
xdat=[]
ydat=[]
zdat=[]
for k in range(0,len(data_dict['0_H_Label'])):
    xdat.append(-(data_dict['0_H_Energy'][k]+0.2))
    ydat.append(data_dict['0_H_Bandgap'][k])
    zdat.append(data_dict['0_H_Dist_max'][k])

p1=plt.scatter(xdat,ydat,c=zdat, vmin=0, vmax=1, s=sdot,marker='o', cmap=cm,edgecolors= "black", label='(H$^+$ +e$^-$) -> H$^*$')
cbar=plt.colorbar(ticks=[0.0, 0.25, 0.5, 0.75, 1.0])
cbar.set_label('Distortion [Ang]', fontsize=size1)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(size2)

xdat2=[]
ydat2=[]
zdat2=[]
for k in range(0,len(data_dict['0_remove_H_Label'])):
    xdat2.append(data_dict['0_remove_H_Energy'][k]+0.2)
    ydat2.append(data_dict['0_remove_H_Bandgap'][k])
    zdat2.append(data_dict['0_remove_H_Dist_max'][k])
p2=plt.scatter(xdat2,ydat2,c=zdat2, vmin=0, vmax=1, s=sdot,marker='s', cmap=cm,edgecolors= "black", label=r'H -> (H$^+$ +e$^-$)')

xdat3=[]
ydat3=[]
zdat3=[]
for k in range(0,len(data_dict['0_OH_Label'])):
    xdat3.append(data_dict['0_OH_Energy'][k]+0.35-0.3)
    ydat3.append(data_dict['0_OH_Bandgap'][k])
    zdat3.append(data_dict['0_OH_Dist_max'][k])
p3=plt.scatter(xdat3,ydat3,c=zdat3, vmin=0, vmax=1, s=sdot,marker='^', cmap=cm,edgecolors= "black", label='H$_2$O -> $^*$OH+(H$^+$ +e$^-$)')

xdat4=[]
ydat4=[]
zdat4=[]
for k in range(0,len(data_dict['0_O_Label'])):
    xdat4.append(data_dict['0_O_Energy'][k]+0.05)
    ydat4.append(data_dict['0_O_Bandgap'][k])
    zdat4.append(data_dict['0_O_Dist_max'][k])
p4=plt.scatter(xdat4,ydat4,c=zdat4, vmin=0, vmax=1, s=sdot,marker='p', cmap=cm,edgecolors= "black",label='H$_2$O -> $^*$O+2(H$^+$ +e$^-$)')


def sigmoid(x,mi, mx): return mi + (mx-mi)*(lambda t: (1+16000**(-t+0.5))**(-8) )( (x-mi)/(mx-mi) )
x = np.linspace(-4,7,1000)
plt.plot(x+0.2-0.6, sigmoid(x, 0.0, 1.5),'b-', lw=3, alpha=0.5)
plt.fill_between(x+0.2-0.6,sigmoid(x, 0.0, 1.5),color='b',alpha=0.2)

def sigmoid(x,mi, mx): return mi + (mx-mi)*(lambda t: (1+16000**(t-0.5))**(-12) )( (x-mi)/(mx-mi) )
x = np.linspace(-4,7,1000)
plt.plot(x-0.3, sigmoid(x, 0.0, 1.2),'r-', lw=3, alpha=0.5)
plt.fill_between(x-0.3,sigmoid(x, 0.0, 1.2),color='r',alpha=0.2)

#plt.legend(fontsize=size2, loc=2, ncol=2)

legend1=plt.legend([p1], [r'(H$^+$ +e$^-$) $\rightarrow$ H$^*$'], loc=2, fontsize=size2)
legend2=plt.legend([p2,p3,p4], [r'H$^*$ $\rightarrow$ (H$^+$ +e$^-$)',r'H$_2$O $\rightarrow$ $^*$OH+(H$^+$ +e$^-$)',r'H$_2$O $\rightarrow$ $^*$O+2(H$^+$ +e$^-$)'], loc=1, fontsize=size2)
plt.gca().add_artist(legend1)
plt.text(-1.8, 2.6, r'$\bf{a}$ MOF-525', ha='left', fontsize=size2+14)

plt.text(-1.48, 0.1, r'Reduction', ha='left', fontsize=size2)
plt.text(1.58, 0.1, r'Oxidation', ha='left', fontsize=size2)


plt.xticks(fontsize=size2)
plt.yticks(fontsize=size2)
plt.xlabel('E / V$_{RHE}$ [eV]', fontsize=size1)
plt.ylabel('Bandgap [eV]', fontsize=size1)
plt.plot([0,0],[0,4],'k-',linewidth=4)
#plt.title(Iron+'Fe with ' +Intermediates)
plt.ylim([0,2.5])
plt.xlim([-1.5,2.0])
#plt.axis([-1, 1, -1, 10])
plt.savefig('0Fe_with_SPIN_potential_map.png', dpi=400, bbox_inches='tight')
plt.show()
     
