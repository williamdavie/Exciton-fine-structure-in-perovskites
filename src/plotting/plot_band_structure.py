
'''
Script to plot the band structure given an input nscf output file.

Band data collected and setup using class: readQEoutput.

William Davie 07/25
'''

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.axes import Axes

import sys,os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))

from QEpostprocessing.read_QE_output import readQEouput
from format import runFormat

#fetch main data

data = readQEouput('CsPbI_0.0_0.0_nscf.out')
kpoints, bandData = data.fetchBandStructure()

# data2 etc .. 

# define high sym points

highsymmetrypoints = {
    
    'M' : np.array([0.5,0.5,0.0]),
    'Γ' : np.array([0,0,0]),
    'X' : np.array([0.5,0.0,0.0]),
    
}

# write labels (for high sym points)

highsymKeys = list(highsymmetrypoints.keys())
label_indices = [0]
labels = [highsymKeys[0]]
for i in range(len(kpoints)):
        
    for j in range(len(highsymKeys)):
            
        a,b,c = highsymmetrypoints[highsymKeys[j]]
        a2,b2,c2 =  kpoints[i]
            
        if (a, b, c) == (a2, b2, c2):
                
            label_indices.append(i)
            labels.append(highsymKeys[j])
            

# used for plotting multiple structures if necessary

def addBandStructure(ax: Axes, kpoints: np.ndarray, bandData: np.ndarray, startband: int=0, endband: int=0,
                     col: str='blue',alpha: float=1,label: str='',linestyle='-'):
    '''
    adds a band structure to an predifined ax.
    '''
    
    numKpoints, numBands = bandData.shape
    
    if endband == 0:
        endband = numBands

    for band in range(startband,endband):
        
        ax.plot(bandData[:,band],c=col,alpha=alpha,label=label,linestyle=linestyle)
    
    return ax

# set up figure:

runFormat()
fig, ax = plt.subplots()
ax.set(ylabel = r'$E(k)$ (eV)')
ax.set_xticks(label_indices)
ax.set_xticklabels(labels)
addBandStructure(ax,kpoints,bandData,startband=119,endband=122)

plt.show()
#plt.savefig('bandstructure.jpg')