import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.axes import Axes

def fetchHighSymLabels(kpoints: np.ndarray, highsymmetrypoints: dict):
    

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
                
    return label_indices, labels
            

# used for plotting multiple structures if necessary

def addBandStructure(ax: Axes, kpoints: np.ndarray, bandData: np.ndarray, efermi: float, startband: int=0, endband: int=0,
                     col: str='blue',alpha: float=1,label: str='',linestyle='-'):
    '''
    adds a band structure to an predifined ax.
    '''
    
    numKpoints, numBands = bandData.shape
    
    if endband == 0:
        endband = numBands

    labeled = False
    for band in range(startband,endband):
        
        if labeled == False:
            ax.plot(bandData[:,band]-efermi,c=col,alpha=alpha,label=label,linestyle=linestyle)
            labeled = True
        else:
            ax.plot(bandData[:,band]-efermi,c=col,alpha=alpha,linestyle=linestyle)
    
    return ax