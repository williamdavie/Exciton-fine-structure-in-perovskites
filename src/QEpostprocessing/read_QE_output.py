'''

Functions to read data from QE .out files

'''

import numpy as np
import math
import re
import os   


class readQEouput():
    
    def __init__(self, filename: str):
    
        self.file = open(filename,'r')
        self.lines = self.file.readlines()
        self.numLines = len(self.lines)
        
    
    def getCelldims(self) -> np.ndarray:
        '''
        Retrieves the input celldims from the output file.
        
        only for ibrav=0
        '''
        
        # written for direct parsing into effective mass calculation

        self.celldims = np.zeros((3))
        scale = 1
        
        for i in range(self.numLines):
            
            line = self.lines[i]
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
            if 'celldm(1)' in line:
        
                scale = float(linesplit[1]) * 0.529177 # Bohr -> Anstong conversion
                
            for i in range(3):
                
                if f'a({i+1})' in line:
                    
                    self.celldims[i] = float(linesplit[3+i]) * scale
                    
        return self.celldims
    

    def fetchBandGap(self) -> float:
        '''
        Fetches the band gap for a QE PWSCF .out file
        '''
        self.gap = 0
        for i in range(self.numLines):
            
            line = self.lines[i]
        
            if "highest occupied, lowest unoccupied level (ev):" in line:
            
                line = re.split('[ ]+',line)
                self.Elumo = float(line[-1])
                self.Ehomo = float(line[-2])
            
                self.gap = abs(self.Ehomo - self.Elumo)
                
                
        return self.gap
    
    
    def fetchBandStructure(self, verbosity: str='low') -> tuple[np.ndarray]:
        '''
        Reads the bandstructure from an nscf output file
        '''
        
        self.nbnds = 0
        self.Nk = 0
        
        assert verbosity == 'high' or verbosity == 'low', 'Incorrect verbosity input'
        
        # first fetch number of bands and k points
        
        for i in range(self.numLines):
            
            line = self.lines[i]
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
            if "number of Kohn-Sham states" in line:
                self.nbnds = int(linesplit[-1])
            
            if "number of k points" in line and verbosity == 'low':
                self.Nk = int(linesplit[-1])
                
            if 'init_us_2' in line and verbosity == 'high': # depends on whether high or low verbosity has been used
                self.Nk = int(linesplit[-2])
                
                
        # then set the values
        self.kpoints = np.zeros((self.Nk,3))
        self.bandData = np.zeros((self.Nk,self.nbnds))
        numBandLines = math.ceil(self.nbnds / 8)
        currentIndex = 0
        
        for i in range(self.numLines):
            
            line = self.lines[i]
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
    
            
            if "bands (ev)" in line:
                
                if verbosity == 'high':
                    
                    if 'occupation numbers' in self.lines[i+3 + numBandLines]:
                        pass # exit out of if statement
                    else:
                        continue # move onto next line
              
                # --------- K points ---------
                
                # minus signs fill space so split doesn't work correctly
                
                minusSigns = [i for i, minusSign in enumerate(line) if minusSign == '-']
                newLine = line
                numFix = 0
                for j in minusSigns:
                    newLine = newLine[:(j+numFix)] + " " + newLine[(j+numFix):]
                    numFix += 1
                
                newLine = newLine.strip()
                newLine = re.split('[ ]+',newLine)
                
                kpoint = np.array([newLine[2],newLine[3],newLine[4]],dtype=np.float32)
            
                self.kpoints[currentIndex] = kpoint
                
                # --------- Energy values ---------
                
                Block = " ".join(line.strip() for line in self.lines[(i+2):(i+2+numBandLines)])
                try:
                    Evals = np.array(re.split('[ ]+',Block))
                except:
                    "file in incorrect format, looking for nscf .out"
                
                try:
                    self.bandData[currentIndex] = Evals
                except:
                    pass

                currentIndex += 1
                
    
        return self.kpoints,self.bandData
    
    
    def fetchHVBandLCB(self) -> tuple[np.ndarray]:
        '''
        Returns the highest valance band and lowest conduction band respectively
        uses results from both functions above: fetchBandstructure and fetchBandgap
        '''    
        
        HVBindex = 0,0
        LCBindex = 0,0

        for i in range(self.Nk):
            for j in range(self.nbnds):
                if self.bandData[i,j] == self.Ehomo:
                    HVBindex = i,j
                    
                if self.bandData[i,j] == self.Elumo:
                    LCBindex = i,j
                    
        #assert HVBindex[0] == LCBindex[0], "Band gap must be direct"        
        if HVBindex[0] != LCBindex[0]:
            print('Warning <!> : Band gap indirect')
        
        HVB = self.bandData[:,HVBindex[1]]
        LCB = self.bandData[:,LCBindex[1]]
        
        return HVB,LCB

        
def fetchDirBandGap(directory: str, 
                    store: bool=True, printResults: bool=True, outputfilename: str='bandGapLandscape'):
    '''
    Calls fetchBandGap() on a set of files in a directory.
    '''
    
    # if store=True function is specifc to the 45 structures we are interested in 
    gapLandscape = np.zeros(shape=(9,9))
    angleVals = np.array([0,2.5,5,7.5,10,12.5,15,17.5,20])
    
    for filename in os.listdir(directory):
        
        fileData = readQEouput(directory + '/' + filename)
        gap = fileData.fetchBandGap()
        
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
            
            
    np.save(outputfilename,gapLandscape)
