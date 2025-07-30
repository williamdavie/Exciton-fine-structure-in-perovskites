'''
Plots the band structure across a kpath Z-Γ-X-M-Γ-X
'''

import matplotlib.pyplot as plt 
import numpy as np
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))
from QEpostprocessing.read_QE_output import readQEouput
from plotting_functions import runFormat,addBandStructure,fetchHighSymLabels,fetchColour


data = readQEouput('./src/Scripts/ScriptData/CsPbI_0.0_0.0_nscf_band.out')
kpoints, bandData = data.fetchBandStructure(verbosity='high')
data_SR = readQEouput('./src/Scripts/ScriptData/CsPbI_0.0_0.0_nscf_SR.out')
kpoints_SR, bandData_SR = data_SR.fetchBandStructure(verbosity='low')

highsymmetrypoints = {
    
    'Z' : np.array([0.0,0.0,0.5]),
    'Γ' : np.array([0,0,0]),
    'X' : np.array([0.5,0.0,0.0]),
    'M' : np.array([0.5,0.5,0.0]),
    'Γ' : np.array([0.0,0.0,0.0]),
    'Y' : np.array([0.0,0.5,0.0]),
    
}
            
label_indicies, labels = fetchHighSymLabels(kpoints,highsymmetrypoints)

runFormat()

label_indicies, labels = fetchHighSymLabels(kpoints,highsymmetrypoints) #  same for all plots 

betastr = r'\beta'

colour = fetchColour('red')

'''
Notes on plotting in report:

Chapter 4:
Fig 1: Band structure + PDOS, figsize = (5,6), xlim = [-10,5]

Fig X: Scalar Rel vs Fully Rel, figsize = (6,4)

'''

fig, ax = plt.subplots(figsize=(5,6)) 

addBandStructure(ax,kpoints_SR,bandData_SR,-4.8629,linewidth=1.5,linestyle='--',col='black',label='SR')
#addBandStructure(ax,kpoints,bandData,-4.6619,linewidth=1.5,cfetchColour('lightred')ol=,label='SOC')


ax.margins(x=0)
ax.set(xlabel=r'$\mathbf{k}$',ylabel = r'$E - E_f$ (eV)')
ax.set_xticks(label_indicies)
ax.set_xticklabels(labels)

print(label_indicies)
print(labels)
ax.set(ylim=[-1,4])
ax.set(xlim=[30,50])

ax.set_title('SR',loc='left')

#ax.legend()
plt.savefig('SRonly.pdf')
plt.show()