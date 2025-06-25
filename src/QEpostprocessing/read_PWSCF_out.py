
'''
A selection of functions to read PWSCF .out files 
'''

import numpy as np
import re
import os    


def fetchBandGap(filename: str) -> float:
    '''
    Attempts to fetch the band gap for a QE .out file
    
    Assumes band gap at Gamma poitn 
    '''
    
    file = open(filename,'r')
    lines = file.readlines()

    numbands = 140
    fermiLevel = 0

    values = np.zeros(numbands)

    for i in range(len(lines)):
        
        if "highest occupied, lowest unoccupied level (ev):" in lines[i]:
            
            line = lines[i].strip()
            line = re.split('[ ]+',line)
            Hop = float(line[-1])
            Lop = float(line[-2])
            
            # guess fermi level as between the two
            
            gap = (Hop - Lop)

    return gap

#------------------------------------------------

def fetchDirBandGap(directory, store=True, printResults=True):
    '''
    Calls fetchBandGap() on a set of files in a directory.
    '''
    
    # if store=True function is specifc to the 45 structures we are interested in 
    
    # sets up bandGap landscape ready for plotting. 
    
    gapLandscape = np.zeros(shape=(9,9))
    angleVals = np.array([0,2.5,5,7.5,10,12.5,15,17.5,20])
    
    for filename in os.listdir(directory):
    
        gap = fetchBandGap(directory + '/' + filename)
        
        if printResults == True:
            print(filename)
            print(gap)
            
        if store == True:
            
            # assumes filename contains _beta_delta_ as setup
            
            underscore = [i for i, underscore in enumerate(filename) if underscore == '_']
            beta = float(filename[underscore[0]+1:underscore[1]])
            delta = float(filename[underscore[1]+1:underscore[2]])
            
            betaIndex = np.where(angleVals==beta)
            deltaIndex = np.where(angleVals==delta)
            
            print(betaIndex)
            print(deltaIndex)
            
            gapLandscape[betaIndex,deltaIndex] = gap
            
    np.save('bandGapLandscape',gapLandscape)
            
        

fetchDirBandGap('./BandGaps_outFiles')