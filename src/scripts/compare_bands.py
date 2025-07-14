'''
This script compares the position of the valance band maximum and conduction band maximum in two cases:

- Compares an undistored structure (0,0) to a beta-only distorted structure (15,0) to test the effect of beta

- Compares a beta-only distorted structure (15,0 )to a `double distored' structure (15,12.5) to test the effect of delta 

William Davie 07/25
'''

# fermi levels were found from scf calculations with occupations smearing:

import matplotlib.pyplot as plt 
import numpy as np
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))
from QEpostprocessing.read_QE_output import readQEouput
from plotting.format import runFormat
from plotting.plot_band_structure import fetchHighSymLabels,addBandStructure

#define and fetch main data

beta = [0,15,15]
delta = [0,0,12.5]
fermiLevels = [-4.1688,-3.6904,-4.2064] 

data = readQEouput('./Scripts/ScriptData/CsPbI_0.0_0.0_nscf.out')
kpoints, bandData = data.fetchBandStructure()
data2 = readQEouput('./Scripts/ScriptData/CsPbI_15.0_0.0_nscf.out')
kpoints2, bandData2 = data2.fetchBandStructure()
data3 = readQEouput('./Scripts/ScriptData/CsPbI_15.0_12.5_nscf.out')
kpoints3, bandData3 = data3.fetchBandStructure()

# define high sym points, see output files 

highsymmetrypoints = {
    
    'M' : np.array([0.5,0.5,0.0]),
    'Γ' : np.array([0,0,0]),
    'X' : np.array([0.5,0.0,0.0]),
    
}
            
label_indicies, labels = fetchHighSymLabels(kpoints,highsymmetrypoints)

runFormat()

label_indicies, labels = fetchHighSymLabels(kpoints,highsymmetrypoints) #  same for all plots 

betastr = r'\beta'

# known index of valance band minimum and conduction maximum, VBM_index = 120,CBM_index = 121

# 1st comparison (0,0 vs 15,0):

fig, ax = plt.subplots()
ax.set(xlabel=r'$\mathbf{k}$',ylabel = r'$E - E_f$ (eV)')
ax.set_xticks(label_indicies)
ax.set_xticklabels(labels)
addBandStructure(ax,kpoints,bandData,fermiLevels[0],startband=119,endband=122,label=f'${betastr} = {beta[0]},\delta = {delta[0]}$',col='#F24405')
addBandStructure(ax,kpoints2,bandData2,fermiLevels[1],startband=119,endband=122,linestyle=':',label=f'${betastr} = {beta[1]},\delta = {delta[1]}$',col='#F24405')
ax.legend()
plt.savefig('beta_vs_undistorted.pdf')

# 2nd comparison ()

fig, ax = plt.subplots()
ax.set(xlabel=r'$\mathbf{k}$',ylabel = r'$E - E_f$ (eV)')
ax.set_xticks(label_indicies)
ax.set_xticklabels(labels)
addBandStructure(ax,kpoints2,bandData2,fermiLevels[1],startband=119,endband=122,label=f'${betastr} = {beta[1]},\delta = {delta[1]}$',col='#F24405')
addBandStructure(ax,kpoints3,bandData3,fermiLevels[2],startband=119,endband=122,linestyle=':',label=f'${betastr} = {beta[2]},\delta = {delta[2]}$',col='#F24405')
ax.legend()
plt.savefig('double_vs_betaonly.pdf')













