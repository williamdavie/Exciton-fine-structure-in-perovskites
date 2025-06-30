
import numpy as np
import math
import re

'''
A selection of functions to read Quantum Espresso output files 
'''

class readQEouput():
    
    def __init__(self, filename: str):
    
        self.file = open(filename,'r')
        self.lines = self.file.readlines()
        self.numLines = len(self.lines)
        
    
    def getCelldims(self):
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
        
        for i in range(self.numLines):
            
            line = self.lines[i]
        
            if "highest occupied, lowest unoccupied level (ev):" in line:
            
                line = re.split('[ ]+',line)
                self.Elumo = float(line[-1])
                self.Ehomo = float(line[-2])
            
                self.gap = abs(self.Ehomo - self.Elumo)
                
                
        return self.gap
    
    
    def fetchBandStructure(self) -> tuple[np.ndarray]:
        '''
        Reads the bandstructure from an nscf output file
        '''
        
        self.nbnds = 0
        self.Nk = 0
        
        # first fetch number of bands and k points
        
        for i in range(self.numLines):
            
            line = self.lines[i]
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
            if "number of Kohn-Sham states" in line:
                self.nbnds = int(linesplit[-1])
            
            if "number of k points" in line:
                self.Nk = int(linesplit[-1])
                
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
                
                self.bandData[currentIndex] = Evals

                currentIndex += 1
                
    
        return self.kpoints,self.bandData
    
    
    def fetchLCBandHVB(self):
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
                    
        assert HVBindex[0] == LCBindex[0], "Band gap must be direct"        

        HVB = self.bandData[:,HVBindex[1]]
        LCB = self.bandData[:,LCBindex[1]]
        
        return HVB,LCB

    
    
    
    
    