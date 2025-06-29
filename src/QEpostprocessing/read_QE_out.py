
import numpy as np
import re

'''
A selection of functions to read Quantum Espresso output files 
'''
#----------------------------------------------------
# Band Gap:
#----------------------------------------------------

def fetchBandGap(filename: str, nbnd: int=140) -> float:
    '''
    Attempts to fetch the band gap for a QE PWSCF .out file
    '''
    
    file = open(filename,'r')
    lines = file.readlines()

    values = np.zeros(nbnd)

    gap = 0 
    for i in range(len(lines)):
        
        if "highest occupied, lowest unoccupied level (ev):" in lines[i]:
            
            line = lines[i].strip()
            line = re.split('[ ]+',line)
            Hop = float(line[-1])
            Lop = float(line[-2])
            
            gap = (Hop - Lop)
    
    file.close()

    return gap

import os    

def fetchDirBandGap(directory, store=True, printResults=True):
    '''
    Calls fetchBandGap() on a set of files in a directory.
    '''
    
    # if store=True function is specifc to the 45 structures we are interested in 
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
            
            
#fetchDirBandGap('./BandGaps_outFiles')

#----------------------------------------------------
# Band Structure:
#----------------------------------------------------

def fetchBandData(filename: str):
    '''
    Reads the data from a Quantum Espresso band structure. 
    '''
    
    file = open(filename,'r')
    lines = file.readlines()
    
    # fetch num bands and num k points
    
    header = lines[0].strip()
    header = re.split('[,=/]',header)
    
    numBands = int(header[1])
    numKpoints = int(header[3])
    
    data = np.zeros((numKpoints,numBands),dtype=np.float32)
    
    print(numBands,numKpoints)
    
    # number of values per line is 10
    
    # calculates the size of divisions. 
    
    if numBands % 10 == 0: 
        div = numBands // 10 + 1 
    else:
        div = numBands // 10 + 2
    
    kvalues = np.zeros((numKpoints,3),dtype=np.float32)
        
    for index, line in enumerate(lines[1:]):
        
        line = line.strip('\n')
        
        currentkpoint = index // div
        divLine = index % div
        
        if divLine == 0: # then k point is here.
            
            line = re.split('[ ]+',line)
            kvalues[currentkpoint]  = np.array([float(line[1]),float(line[2]),float(line[3])])
            
        else:
            
            line = re.split('[ ]+',line)
            
            # computes the start and end points of the current line
            a = (divLine - 1) * 10
            b = a + len(line[1:11])
            
            data[currentkpoint][a:b] = line[1:11]
        
    file.close()
    
    return kvalues, data
    
   