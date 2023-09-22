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
from sklearn.metrics import mean_absolute_error as mae

# From Fend
#ESlab=-1908.364500 # Need PW500 number

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
metals=['Co', 'Pd', 'Ni', 'Pt', 'Ir', 'Cu','Zn','Mn','Rh','Fe']

plt.figure()
x=[]
y=[]
Intermed='H'
for Intermed in intermediates: 
    plt.figure()
    plt.locator_params(axis='x',nbins=5);plt.locator_params(axis='y',nbins=5);
    x=[]
    y=[]
    for metal in metals:
            index_x=MOF525_dict[Intermed+'_Metal'].index(metal)
            index_y=porphyrin_slab_dict[Intermed+'_Metal'].index(metal)
            plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_slab_dict[Intermed+'_Energy'][index_y],color='b')
            plt.text(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_slab_dict[Intermed+'_Energy'][index_y],porphyrin_slab_dict[Intermed+'_Metal'][index_y])
            
            x.append(MOF525_dict[Intermed+'_Energy'][index_x])
            y.append(porphyrin_slab_dict[Intermed+'_Energy'][index_y])
    
    plt.annotate('slab vs MOF525: '+Intermed, xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                         bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                         horizontalalignment='left', verticalalignment='bottom')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn=np.min(x)
    mx=np.max(x)
    x1=np.linspace(mn,mx,500)
    y1=gradient*x1+intercept
    y2=gradient*x1+intercept
    plt.plot(x1,y1,'-b')
    plt.text(x1[-1]-1.5,y1[-1]-0.5,str(np.round(gradient,1))+'x+'+str(np.round(intercept,1))+'\n $r^2$=' +str(np.round(r_value,1)))
                
    plt.xticks(fontsize=size2)
    plt.yticks(fontsize=size2)
    plt.xlabel('MOF525',fontsize=size1)
    plt.ylabel('Porphyrin slab',fontsize=size1)
    #plt.savefig('MOF525_Metals_COOH_vs_H.png', dpi=400, bbox_inches='tight')
    

plt.figure()
x=[]
y=[]
Intermed='H'
for Intermed in intermediates: 
    plt.figure()
    plt.locator_params(axis='x',nbins=5);plt.locator_params(axis='y',nbins=5);
    x=[]
    y=[]
    for metal in metals:
            index_x=MOF525_dict[Intermed+'_Metal'].index(metal)
            index_y=porphyrin_mol_dict[Intermed+'_Metal'].index(metal)
            plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_mol_dict[Intermed+'_Energy'][index_y],color='b')
            plt.text(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_mol_dict[Intermed+'_Energy'][index_y],porphyrin_mol_dict[Intermed+'_Metal'][index_y])
            
            x.append(MOF525_dict[Intermed+'_Energy'][index_x])
            y.append(porphyrin_mol_dict[Intermed+'_Energy'][index_y])
    
    plt.annotate('Mol vs MOF525: '+Intermed, xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                         bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                         horizontalalignment='left', verticalalignment='bottom')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn=np.min(x)
    mx=np.max(x)
    x1=np.linspace(mn,mx,500)
    y1=gradient*x1+intercept
    plt.plot(x1,y1,'-b')
    plt.text(x1[-1]-1.5,y1[-1]-0.5,str(np.round(gradient,1))+'x+'+str(np.round(intercept,1))+'\n $r^2$=' +str(np.round(r_value,1)))
                
    plt.xticks(fontsize=size2)
    plt.yticks(fontsize=size2)
    plt.xlabel('MOF-525',fontsize=size1)
    plt.ylabel('Porphyrin mol',fontsize=size1)
    #plt.savefig('MOF525_Metals_COOH_vs_H.png', dpi=400, bbox_inches='tight')


#corr=[[0,0],[1,0.5,0.25]]
#sns.heatmap(corr, vmin=-1.0, vmax=1.0, square=True, cmap="RdBu")

m=[]
for metal in metals:
    m.append(metal)
    
d = {'Metals': m}
dfm = pd.DataFrame(data=d)

plt.figure()
x=[]
y=[]
z=[]
Intermed='H'
for Intermed in ['H', 'COOH', 'OH', 'OOH','O']: 
    plt.figure()

    if Intermed=='H':
        plt.annotate(r'$\bf{a}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
    elif Intermed=='COOH':
        plt.annotate(r'$\bf{b}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')  

    elif Intermed=='OH':
        plt.annotate(r'$\bf{c}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom') 

    elif Intermed=='OOH':
        plt.annotate(r'$\bf{d}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom') 

    elif Intermed=='O':
        plt.annotate(r'$\bf{e}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom') 

    plt.locator_params(axis='x',nbins=5);plt.locator_params(axis='y',nbins=5);
    x=[]
    y=[]
    z=[]
    m=[]
    for metal in metals:
            index_x=MOF525_dict[Intermed+'_Metal'].index(metal)
            index_y=porphyrin_mol_dict[Intermed+'_Metal'].index(metal)
            index_z=porphyrin_slab_dict[Intermed+'_Metal'].index(metal)
            plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_mol_dict[Intermed+'_Energy'][index_y],color='b')
            plt.text(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_mol_dict[Intermed+'_Energy'][index_y],porphyrin_mol_dict[Intermed+'_Metal'][index_y],fontsize=size2-4)

            plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_slab_dict[Intermed+'_Energy'][index_z],color='r')
            #plt.text(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_slab_dict[Intermed+'_Energy'][index_z],porphyrin_slab_dict[Intermed+'_Metal'][index_z])

            x.append(MOF525_dict[Intermed+'_Energy'][index_x])
            y.append(porphyrin_mol_dict[Intermed+'_Energy'][index_y])
            z.append(porphyrin_slab_dict[Intermed+'_Energy'][index_z])
            m.append(metal)
    
    d = {'Metals': m, 'MOF-525_'+Intermed: x, 'Mol_'+Intermed: y, 'Slab_'+Intermed: z}
    df = pd.DataFrame(data=d)
    dfm=pd.merge(dfm, df, on='Metals', how='outer')
    plt.annotate(Intermed+' vs. '+Intermed, xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                         bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                         horizontalalignment='left', verticalalignment='bottom')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    mn=np.min(x)
    mx=np.max(x)
    x1=np.linspace(mn,mx,500)
    y1=gradient*x1+intercept
    plt.plot(x1,y1,'-b')
    #plt.text(mn+0.0,np.max(y1)-0.2*(np.max(y1)-np.min(y1)),r'TCCP$_{mol}$'+'\n' + str(np.round(gradient,1))+'MOF525+'+str(np.round(intercept,1))+'\n'+'$r^2$=' +str(np.round(r_value,1)),c='b', fontsize=size2-4)
         
    y2=gradient*np.asarray(x)+intercept
    error = mae(y, y2)
    plt.annotate(r'TCCP$_{mol}$'+'\n' + str(np.round(gradient,1))+'MOF-525+'+str(np.round(intercept,1))+'\n'+'$r^2$=' +str(np.round(r_value,1))+', MAE=' +str(np.round(error,2)), xy=(0.03, 0.6), xycoords='axes fraction', fontsize=size2-4,c='b',horizontalalignment='left', verticalalignment='bottom')

    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,z)
    mn=np.min(x)
    mx=np.max(x)
    x1=np.linspace(mn,mx,500)
    y1=gradient*x1+intercept
    plt.plot(x1,y1,'-r')
#    plt.text(mn+1.0,np.min(y1)+0.2,r'TCCP$_{slab}$'+'\n'+str(np.round(gradient,1))+'MOF525+'+str(np.round(intercept,1))+'\n'+'$r^2$=' +str(np.round(r_value,1)),c='r',fontsize=size2-4)
    y2=gradient*np.asarray(x)+intercept
    error = mae(z, y2)
    plt.annotate(r'TCCP$_{slab}$'+'\n'+str(np.round(gradient,1))+'MOF-525+'+str(np.round(intercept,1))+'\n'+'$r^2$=' +str(np.round(r_value,1))+', MAE=' +str(np.round(error,2)), xy=(0.53, 0.1), xycoords='axes fraction', fontsize=size2-4,c='r',horizontalalignment='left', verticalalignment='bottom')
    

    plt.xticks(fontsize=size2)
    plt.yticks(fontsize=size2)
    plt.xlabel('$\Delta$E MOF-525',fontsize=size2)
    plt.ylabel('$\Delta$E{'+'TCCP$_{mol}$'+', '+'TCCP$_{slab}$}',fontsize=size2)
    plt.savefig('MOF525_Metals_'+Intermed+'_'+Intermed+'.png', dpi=400, bbox_inches='tight')

plt.figure()
import seaborn as sns
corr = dfm.select_dtypes('number').corr()
# plot the heatmap
#sns.heatmap(corr)
sns.heatmap(corr, vmin=0.7, vmax=1.0, square=True, cmap="RdBu")



# Volcano H 
top=-0.2    

plt.figure()
plt.plot([-1,top],[-1-top,0],'k',lw=2,zorder=1)
plt.plot([top,1],[0,-1+top],'k',lw=2,zorder=1)
plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=4)
s1=50
size3=size2-10
Intermed='H'
plt.annotate(r'$\bf{a}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
for metal in metals:
    index_x=MOF525_dict[Intermed+'_Metal'].index(metal)
    index_y=porphyrin_mol_dict[Intermed+'_Metal'].index(metal)
    index_z=porphyrin_slab_dict[Intermed+'_Metal'].index(metal)
    
    if MOF525_dict[Intermed+'_Energy'][index_x] > top and MOF525_dict[Intermed+'_Energy'][index_x] < 0.7:
        p1=plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],-(MOF525_dict[Intermed+'_Energy'][index_x]-top),color='k', marker='s', s=s1)
        plt.text(MOF525_dict[Intermed+'_Energy'][index_x]+0.03,-(MOF525_dict[Intermed+'_Energy'][index_x]-top),metal,fontsize=size3)
    elif MOF525_dict[Intermed+'_Energy'][index_x] < top:
        plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],MOF525_dict[Intermed+'_Energy'][index_x]-top,color='k', marker='s', s=s1)  
        plt.text(MOF525_dict[Intermed+'_Energy'][index_x]-0.08,(MOF525_dict[Intermed+'_Energy'][index_x]-top),metal,fontsize=size3)

    if porphyrin_mol_dict[Intermed+'_Energy'][index_y] > top and porphyrin_mol_dict[Intermed+'_Energy'][index_y] < 0.7:
        p2=plt.scatter(porphyrin_mol_dict[Intermed+'_Energy'][index_y],-(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-top),color='b',marker='o',s=s1)
        plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-0.1,-(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-top)-0.05,metal, color='b',fontsize=size3)
    elif porphyrin_mol_dict[Intermed+'_Energy'][index_y] < top:
        plt.scatter(porphyrin_mol_dict[Intermed+'_Energy'][index_y],porphyrin_mol_dict[Intermed+'_Energy'][index_y]-top,color='b', marker='o',s=s1)
        plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.05,(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-top),metal, color='b',fontsize=size3)
    
    
    if porphyrin_slab_dict[Intermed+'_Energy'][index_z] > top and porphyrin_slab_dict[Intermed+'_Energy'][index_z] < 0.7: 
        p3=plt.scatter(porphyrin_slab_dict[Intermed+'_Energy'][index_z],-(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-top),color='r', marker='v',s=s1)
        plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.01,-(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-top)+0.03,metal, color='r',fontsize=size3)
    elif porphyrin_slab_dict[Intermed+'_Energy'][index_z] < top:
        plt.scatter(porphyrin_slab_dict[Intermed+'_Energy'][index_z],porphyrin_slab_dict[Intermed+'_Energy'][index_z]-top,color='r', marker='v',s=s1)
        plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-0.08,(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-top),metal, color='r',fontsize=size3)
        
        
        
    #plt.text(MOF525_dict[Intermed+'_Energy'][index_x],porphyrin_slab_dict[Intermed+'_Energy'][index_z],porphyrin_slab_dict[Intermed+'_Metal'][index_z])
            
    #plt.scatter(xdata,np.asarray(xdata)-top,c=xcolor, vmin=0, vmax=0.5, s=sdot, cmap=cm,edgecolors= "black",zorder=2)


#plt.annotate("Fe", xy=(data_dict[Iron+'_H_Fe_Energy'],data_dict[Iron+'_H_Fe_Energy']-top), xytext=(data_dict[Iron+'_H_Fe_Energy'],data_dict[Iron+'_H_Fe_Energy']-top-0.5),fontsize=size2-4,arrowprops=dict(arrowstyle="->"))

plt.annotate('HER', xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                         bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                         horizontalalignment='left', verticalalignment='bottom')

legend1=plt.legend([p1,p2,p3], [r'MOF-525',r'TCCP$_{mol}$', r'TCCP$_{slab}$'], loc=1, fontsize=size2-8)
plt.xticks(fontsize=size2)
plt.yticks(fontsize=size2)
plt.xlabel('$\Delta$E$_{H^*}$ [eV]',fontsize=size1)
plt.ylabel('$\eta$ [eV]',fontsize=size1)
#plt.title(Iron+'Fe with ' +Intermediates)
#plt.ylim([0,1.5])
#plt.axis([-1, 1, -1, 10])

plt.xlim([-1.0,1.0])
plt.ylim([-1.0,0.1])
#plt.savefig('Plot_H_added.png', dpi=600, bbox_inches='tight')
plt.savefig('HER_volcano.png',dpi=900, bbox_inches='tight')



#ORR
top=(4.92-3.2)/2.0

plt.figure()
plt.plot([-2,(4.92-3.2)/2.0],[-2,(4.92-3.2)/2.0],c='k',lw=2,zorder=1)
plt.plot([(4.92-3.2)/2.0,3],[4.92-3.2-(4.92-3.2)/2.0,4.92-3.2-3],c='k',lw=2,zorder=1)
plt.plot([-2,3],[1.23,1.23],c='b',lw=2)
plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
    

Intermed='OH'
plt.annotate(r'$\bf{b}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
for metal in metals:
    index_x=MOF525_dict[Intermed+'_Metal'].index(metal)
    index_y=porphyrin_mol_dict[Intermed+'_Metal'].index(metal)
    index_z=porphyrin_slab_dict[Intermed+'_Metal'].index(metal)
    
    if MOF525_dict[Intermed+'_Energy'][index_x] +0.35-0.3> top:
        p1=plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x]+0.35-0.3,2*top-(MOF525_dict[Intermed+'_Energy'][index_x]+0.35-0.3),color='k', marker='s', s=s1)
        plt.text(MOF525_dict[Intermed+'_Energy'][index_x]+0.03,2*top-(MOF525_dict[Intermed+'_Energy'][index_x]+0.35-0.3)+0.2,metal,fontsize=size3)
    elif MOF525_dict[Intermed+'_Energy'][index_x] +0.35-0.3 < top:
        plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x]+0.35-0.3,MOF525_dict[Intermed+'_Energy'][index_x]+0.35-0.3,color='k', marker='s', s=s1)  
        plt.text(MOF525_dict[Intermed+'_Energy'][index_x]-0.05,(MOF525_dict[Intermed+'_Energy'][index_x]+0.35-0.3)+0.2,metal,fontsize=size3)
        
    if porphyrin_mol_dict[Intermed+'_Energy'][index_y] +0.35-0.3> top:
        p2=plt.scatter(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.35-0.3,2*top-(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.35-0.3),color='b', marker='o', s=s1)
        plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-0.08,2*top-(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.35-0.3)-0.13,metal,color='b',fontsize=size3)
    elif porphyrin_mol_dict[Intermed+'_Energy'][index_y] +0.35-0.3 < top:
        plt.scatter(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.35-0.3,porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.35-0.3,color='b', marker='o', s=s1)  
        plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.08,(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.35-0.3)-0.13,metal,color='b',fontsize=size3)

        
    if porphyrin_slab_dict[Intermed+'_Energy'][index_z] +0.35-0.3> top:
        p3=plt.scatter(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.35-0.3,2*top-(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.35-0.3),color='r', marker='v', s=s1)
        plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.03,2*top-(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.35-0.3)+0.05,metal,color='r',fontsize=size3)
    elif porphyrin_slab_dict[Intermed+'_Energy'][index_z] +0.35-0.3 < top:
        plt.scatter(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.35-0.3,porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.35-0.3,color='r', marker='v', s=s1)  
        plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-0.03,(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.35-0.3)+0.05,metal,color='r',fontsize=size3)


    #plt.annotate("Fe", xy=(data_dict[Iron+'_H_Fe_Energy'],data_dict[Iron+'_H_Fe_Energy']-top), xytext=(data_dict[Iron+'_H_Fe_Energy'],data_dict[Iron+'_H_Fe_Energy']-top-0.5),fontsize=size2-4,arrowprops=dict(arrowstyle="->"))

plt.annotate('ORR', xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                         bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                         horizontalalignment='left', verticalalignment='bottom')

legend1=plt.legend([p1,p2,p3], [r'MOF-525',r'TCCP$_{mol}$', r'TCCP$_{slab}$'], loc=1, fontsize=size2-8)
plt.xticks([0,1,2,3],['0.0','1.0','2.0','3.0'],fontsize=size2)
plt.yticks([-1,0,1],['-1.0','0.0','1.0'],fontsize=size2)
plt.xlabel('$\Delta$G$_{~^*\!OH}$ [eV]',fontsize=size1)
plt.ylabel('$\eta$ [eV]',fontsize=size1)
    #plt.title(Iron+'Fe with ' +Intermediates)
    #plt.ylim([0,1.5])
    #plt.axis([-1, 1, -1, 10])

plt.xlim([0.0,3.0])
plt.ylim([-1.3,1.4])

    #plt.savefig('ORR_'+Iron+'Fe_with_BG_SPIN.png', dpi=400, bbox_inches='tight')
plt.savefig('ORR_volcano.png',dpi=900, bbox_inches='tight')


####### writing faradaic efficiencie
dictM={}
### H2 Gas ####

folder='../../../../PhD_project/Formaldehyde_reduction/plotting/'

with open(folder+'data/H2_BEEF-vdw_energy.txt', 'r') as f:
    EH2mBEEF = [line.rstrip('\n') for line in f]
EH2BEEF = [float(i) for i in EH2mBEEF]

# load in slab energies
for metal in ['Cu','Ag','Au','Ir','Pt','Pd','Rh','Sn','In','Pb','Cd','Zn','Ga','Tl','Fe','Ni','Ti','Hg','Ru']:
    with open(folder+'data/'+metal+'_larger_clean.txt', 'r') as f:
        EH2mBEEF = [line.rstrip('\n') for line in f]
    data = [float(i) for i in EH2mBEEF]
    dictM[metal+'_slab']=np.asarray(data)
    
for metal in ['Cu','Ag','Au','Ir','Pt','Pd','Rh','Sn','In','Pb','Cd','Zn','Ga','Tl','Fe','Ni','Ti','Hg','Ru']: 
    with open(folder+'data/'+metal+'_H_larger_clean.txt', 'r') as f:
        EH2mBEEF = [line.rstrip('\n') for line in f]
    data = [float(i) for i in EH2mBEEF]
    dictM[metal+'_slab_H']=np.asarray(data)
    

metal='Au'
val_Au_H=(dictM[metal+'_slab_H']-dictM[metal+'_slab']-0.5*np.asarray(EH2BEEF)).mean()

dictM['Pb_HER']=5.0 ; dictM['Pb_total']=102.4
dictM['Hg_HER']=0.0 ; dictM['Hg_total']=99.5
dictM['Tl_HER']=6.2 ; dictM['Tl_total']=101.3
dictM['In_HER']=3.3 ; dictM['In_total']=100.3
dictM['Sn_HER']=4.6 ; dictM['Sn_total']=100.1
dictM['Cd_HER']=9.4 ; dictM['Cd_total']=103.0
dictM['Hg_HER']=0.0 ; dictM['Hg_total']=99.5

dictM['Au_HER']=10.2 ; dictM['Au_total']=98.0
dictM['Ag_HER']=12.4 ; dictM['Ag_total']=94.6
dictM['Zn_HER']=9.9 ; dictM['Zn_total']=95.4
dictM['Pd_HER']=26.2 ; dictM['Pd_total']=60.2
dictM['Cu_HER']=20.5 ; dictM['Cu_total']=103.5

dictM['Ni_HER']=88.9 ; dictM['Ni_total']=92.4
dictM['Fe_HER']=94.8 ; dictM['Fe_total']=94.8
dictM['Pt_HER']=95.7 ; dictM['Pt_total']=95.8
dictM['Ti_HER']=99.7 ; dictM['Ti_total']=99.7
dictM['Ga_HER']=79.0 ; dictM['Ga_total']=102.0

A1=0.5;

fig=plt.figure()
plt.annotate(r'$\bf{c}$', xy=(-0.1, 1.01), xycoords='axes fraction', fontsize=size1,c='k',horizontalalignment='left', verticalalignment='bottom')
plt.locator_params(axis='x',nbins=5);plt.grid(True);plt.locator_params(axis='y',nbins=5)
plt.xticks(fontsize = size2); plt.yticks(fontsize = size2)
ax=fig.gca()
MetalmeanH_COOH=[]; MetalmeanFara=[]; PorpyrinmeanH=[]

for metal in ['Pd','Pt','Cd','Sn','In','Pb','Cu','Ag','Zn','Tl','Fe','Ni','Hg','Au']:
    if metal == 'Ir' or metal=='Pt' or metal=='Pd' or metal=='Rh' or metal=='Ni' or metal=='Fe':
        GiveColor='r'
    elif metal=='Sn' or metal=='In' or metal=='Pb' or metal=='Cd' or metal=='Tl' or metal=='Hg':
        GiveColor='y'
    elif metal=='Ag' or metal=='Au' or metal=='Zn':
        GiveColor='b'
    elif metal=='Cu':
        GiveColor='c'

    Exval1=dictM[metal+'_slab_H']-dictM[metal+'_slab']-0.5*np.asarray(EH2BEEF)
    Eyval=dictM[metal+'_total']-dictM[metal+'_HER']
    Exval=Exval1
    MetalmeanH_COOH.append(Exval.mean())  
    MetalmeanFara.append(Eyval)
    p4=plt.scatter(np.asarray(Exval).mean(), np.asarray(Eyval).mean(), c=GiveColor,s=300,alpha=0.5)
    if metal =='Cd':
         ax.annotate(metal, xy=(np.asarray(Exval).mean(),np.asarray(Eyval).mean()), xytext=(np.asarray(Exval).mean()+0.1,np.asarray(Eyval).mean()-10),color='black',alpha=A1, fontsize=size2, arrowprops=dict(arrowstyle="-"))

    elif metal=='Pt':
         ax.annotate(metal, xy=(np.asarray(Exval).mean(),np.asarray(Eyval).mean()), xytext=(np.asarray(Exval).mean()+0.3,np.asarray(Eyval).mean()+8),color='black',alpha=A1, fontsize=size2, arrowprops=dict(arrowstyle="-"))
    elif metal =='Ag':
         ax.annotate(metal, xy=(np.asarray(Exval).mean(),np.asarray(Eyval).mean()), xytext=(np.asarray(Exval).mean()-0.05,np.asarray(Eyval).mean()-12),color='black',alpha=A1, fontsize=size2, arrowprops=dict(arrowstyle="-"))
    elif metal =='Cu':
         ax.annotate(metal, xy=(np.asarray(Exval).mean(),np.asarray(Eyval).mean()), xytext=(np.asarray(Exval).mean()+0.03,np.asarray(Eyval).mean()-12),color='black',alpha=A1, fontsize=size2, arrowprops=dict(arrowstyle="-"))
   

plt.xlim([-0.5,1.5])
plt.ylim([0.0,100])

from scipy.optimize import curve_fit
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))
    return (y)
p0 = [max(MetalmeanFara), np.median(MetalmeanH_COOH),1,min(MetalmeanFara)] # this is an mandatory initial guess
popt, pcov = curve_fit(sigmoid, MetalmeanH_COOH, MetalmeanFara,p0)
x = np.linspace(np.min(MetalmeanH_COOH), np.max(MetalmeanH_COOH), 1000)
y = sigmoid(x, *popt)
plt.plot(x,y,'r-',linewidth=4, label='Sigmoid Fit',zorder=0,alpha=0.2)
ax.annotate('Fit to metals', xy=(-0.1,50), xytext=(0.0,50.0),color='red',alpha=A1, fontsize=size2-8, arrowprops=dict(arrowstyle="-", color='red'))


MetalmeanH_COOH=[]; MetalmeanFara=[]; PorpyrinmeanH=[]
    
for metal in ['Sn','In','Pb','Cd','Au','Ag','Zn','Cu','Tl','Hg']:
    if metal == 'Ir' or metal=='Pt' or metal=='Pd' or metal=='Rh':
        GiveColor='r'
    elif metal=='Sn' or metal=='In' or metal=='Pb' or metal=='Cd' or metal=='Tl' or metal=='Hg':
        GiveColor='y'
    elif metal=='Ag' or metal=='Au' or metal=='Zn':
        GiveColor='b'
    elif metal=='Cu':
        GiveColor='c'

    #Ezval1 =dictM[metal+'_COOH']-dictM[metal+'_slab']-np.asarray(ECO2)-0.5*np.asarray(EH2BEEF)
    Exval1=dictM[metal+'_slab_H']-dictM[metal+'_slab']-0.5*np.asarray(EH2BEEF)
    Eyval=dictM[metal+'_total']-dictM[metal+'_HER']
    Exval=Exval1
    MetalmeanH_COOH.append(Exval.mean())  
    MetalmeanFara.append(Eyval)

fit = np.polyfit(MetalmeanH_COOH,MetalmeanFara,1)
fit_fn = np.poly1d(fit) 
plt.plot([-1.0, 2.5],[fit_fn(-1.0),fit_fn(2.5)],'b--',linewidth=4,label='Linear Fit',zorder=0, alpha=0.2)

plt.legend(loc=4, fontsize=size1)


Intermed='H'
for metal in metals:
    index_x=MOF525_dict[Intermed+'_Metal'].index(metal)
    index_y=porphyrin_mol_dict[Intermed+'_Metal'].index(metal)
    index_z=porphyrin_slab_dict[Intermed+'_Metal'].index(metal)
    

    if MOF525_dict[Intermed+'_Energy'][index_x] < 1.5 and MOF525_dict[Intermed+'_Energy'][index_x] > -0.5:
        p1=plt.scatter(MOF525_dict[Intermed+'_Energy'][index_x],sigmoid(MOF525_dict[Intermed+'_Energy'][index_x], *popt),color='k', marker='s', s=s1,zorder=1)
        
        if metal=='Mn':
            plt.text(MOF525_dict[Intermed+'_Energy'][index_x]+0.03,sigmoid(MOF525_dict[Intermed+'_Energy'][index_x], *popt)+0.5,metal,fontsize=size3)
        elif metal=='Ni':
            plt.text(MOF525_dict[Intermed+'_Energy'][index_x]-0.07,sigmoid(MOF525_dict[Intermed+'_Energy'][index_x], *popt)+0.5,metal,fontsize=size3)
        else:
            plt.text(MOF525_dict[Intermed+'_Energy'][index_x]+0.03,sigmoid(MOF525_dict[Intermed+'_Energy'][index_x], *popt),metal,fontsize=size3)
   
    if porphyrin_mol_dict[Intermed+'_Energy'][index_y] < 1.5 and porphyrin_mol_dict[Intermed+'_Energy'][index_y] > -0.5:
        p2=plt.scatter(porphyrin_mol_dict[Intermed+'_Energy'][index_y],sigmoid(porphyrin_mol_dict[Intermed+'_Energy'][index_y], *popt),color='b', marker='o', s=s1,zorder=1)
         
        if metal=='Co':
             plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-0.12,sigmoid(porphyrin_mol_dict[Intermed+'_Energy'][index_y], *popt)-1,metal,fontsize=size3)
        elif metal=='Mn' or metal=='Ni':
             plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-0.03,sigmoid(porphyrin_mol_dict[Intermed+'_Energy'][index_y], *popt)-5,metal,fontsize=size3)
        elif metal=='Fe':
             plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y],sigmoid(porphyrin_mol_dict[Intermed+'_Energy'][index_y], *popt)+2,metal,fontsize=size3)
        elif metal=='Pt':
             plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]-0.07,sigmoid(porphyrin_mol_dict[Intermed+'_Energy'][index_y], *popt)+2,metal,fontsize=size3)
        else:
             plt.text(porphyrin_mol_dict[Intermed+'_Energy'][index_y]+0.03,sigmoid(porphyrin_mol_dict[Intermed+'_Energy'][index_y], *popt),metal,fontsize=size3)
    
    if porphyrin_slab_dict[Intermed+'_Energy'][index_z] < 1.5 and porphyrin_slab_dict[Intermed+'_Energy'][index_z] > -0.5:
        p3=plt.scatter(porphyrin_slab_dict[Intermed+'_Energy'][index_z],sigmoid(porphyrin_slab_dict[Intermed+'_Energy'][index_z], *popt),color='r', marker='v', s=s1,zorder=1)
        
        if metal=='Co':
            plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.05,sigmoid(porphyrin_slab_dict[Intermed+'_Energy'][index_z], *popt)-1,metal,fontsize=size3)
        elif metal=='Fe' or metal=='Ni':
            plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-0.07,sigmoid(porphyrin_slab_dict[Intermed+'_Energy'][index_z], *popt)+2,metal,fontsize=size3)
        elif metal=='Mn':
            plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]-0.07,sigmoid(porphyrin_slab_dict[Intermed+'_Energy'][index_z], *popt)+2,metal,fontsize=size3)
        else:
            plt.text(porphyrin_slab_dict[Intermed+'_Energy'][index_z]+0.03,sigmoid(porphyrin_slab_dict[Intermed+'_Energy'][index_z], *popt),metal,fontsize=size3)


plt.annotate('CO2RR', xy=(0.03, 0.9), xycoords='axes fraction', fontsize=size2-4,
                         bbox=dict(boxstyle='round',facecolor='white', alpha=1.0),
                         horizontalalignment='left', verticalalignment='bottom')

legend1=plt.legend([p1,p2,p3,p4], [r'MOF-525',r'TCCP$_{mol}$', r'TCCP$_{slab}$','Metals'], loc=4, fontsize=size2-8)
plt.xticks(fontsize=size2)
plt.yticks(fontsize=size2)


plt.xlabel('$\Delta$E$_{H^*}$  [eV]',fontsize=size1)
plt.ylabel('CO2RR FE [$\%$]',fontsize=size1)

plt.savefig('CO2RR_FE.png',dpi=900, bbox_inches='tight')



