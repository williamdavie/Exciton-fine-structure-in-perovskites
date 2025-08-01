

'''
Constructs a `double model' unit cell for a distorted 2D Hallide pervoskite with a caesium cation.

General formula: (Cs)_2 (B X_4)

Note that this cell is build via emperical approximations and should be relaxed via DFT. 

William Davie 17/04/25
'''

import numpy as np
from typing import List
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))
from perovskite_2D import perovskite2D

#------------------------------------------------

class BuildUnitCell(perovskite2D):
    
    def __init__(self, celldims: np.ndarray, beta: float, delta: float, outdir: str,
                 centred: bool=True):
        '''
        celldims: (a,b,c)
        beta: distortion angle 1 
        delta: distortion angle 2
        outdir: output directory for files generated 
        centred: is the monolayer in the centre of the unit cell?
        '''
        
        super().__init__(B='Pb', Bmass=207.2, X='I', Xmass='126.9') # change accordingly if working with a different structure.
        
        self.beta = beta*np.pi/180
        self.delta = delta*np.pi/180
        
        self.celldims = celldims
        
        self.outputfile = outdir + f'/Cs{self.B}{self.X}_{beta}_{delta}'
        
        self.espressoOutdir = '/local/data/public/wd324/QEout/' # Change when needed
        
        self.XAtomList = []
        self.BAtomList = []
        self.CsAtomList = []
        
        if centred:
            self.Z = 0.5
        else:
            self.Z = 0.0
        

    def locateBXAtoms(self) -> None:

        '''
        Defines the position of B and X atoms: 

            (X6)      (X3) (X8)
            /          |  /
        B1 - (X1) - B2 - (X2)   
        /          / |
        (X5)     (X7)  (X4)
    
        '''
        
        assert self.beta >= self.delta, 'Beta cannot be smaller than delta.'
        
        # Positions of B-site cation are fixed.
        
        B1pos = np.array([0.0,0.0,self.Z])*self.celldims
        B2pos = np.array([0.5,0.5,self.Z])*self.celldims
        
        # in plane X atom calculation:
        
        # initial positions:
        
        X1pos = np.array([0.25,0.25,self.Z])*self.celldims
        X2pos = np.array([0.75,0.75,self.Z])*self.celldims
        X3pos = np.array([0.25,0.75,self.Z])*self.celldims
        X4pos = np.array([0.75,0.25,self.Z])*self.celldims

        # B - B bond length 

        B_BLength = np.sqrt(self.celldims[0]**2 + self.celldims[1]**2)
        
        cosTheta = (self.celldims[1]/2) / (B_BLength/2) # could also use celldims[0] but easier option is to just switch 0 and 1 within input


        # based on the constraints of delta and beta on the Pb - I - Pb bond, we arrive at:
        
        L = np.sqrt( ( ((1/4 * B_BLength)**2)*(1 + np.tan(self.beta)**2) ) / ( cosTheta**2 * np.tan(self.delta)**2 + 1 ))
        print(L)
        # using this we can easily compute changes to init positions

        dz = cosTheta * L * np.tan(self.delta)
        
        if self.beta != 0:
            dxy = np.sqrt((1/4 * B_BLength * np.tan(self.beta))**2 - dz**2)
        else:
            dxy = 0

        # need two change vectors when cell is rectangular
        
        changeVector1 = dxy * 1/B_BLength * np.array([self.celldims[1],self.celldims[0]])
        changeVector2 = dxy * 1/B_BLength * np.array([-self.celldims[1],self.celldims[0]])
        
        dx1 = changeVector1[0]
        
        dy1 = changeVector1[1]
        
        dx2 = changeVector2[0]
        
        dy2 = changeVector2[1]

        # induce the distortion. 
        
        X1pos += np.array([-dx1,dy1,-dz])
        X2pos += np.array([dx1,-dy1,dz])
        X3pos += np.array([-dx2,dy2,dz])
        X4pos += np.array([dx2,-dy2,-dz])
        
        # out of plane calculation:
        # assuming perpendiclar to plane generate by in plane atoms
        
        outVector = np.cross(B2pos - X1pos,B2pos - X4pos)

        outUnitVector = 1/np.linalg.norm(outVector) * outVector
        
        print(L/np.cos(self.delta))
        
        BondL = np.sqrt( ( (1/4 * B_BLength**2)*(1 + np.tan(self.beta)**2) ) ) / 2


        X5pos = np.array([outUnitVector[0],-outUnitVector[1],outUnitVector[2]]) * BondL + B1pos
        X6pos = np.array([-outUnitVector[0],outUnitVector[1],-outUnitVector[2]]) * BondL + B1pos
        X7pos = outUnitVector * BondL + B2pos
        X8pos = np.array([-outUnitVector[0],-outUnitVector[1],-outUnitVector[2]])  * BondL + B2pos
        
        # rescale
        
        self.BAtomList = [B/self.celldims for B in [B1pos,B2pos]]
        
        XAtomList = [X1pos,X2pos,X3pos,X4pos,X5pos,X6pos,X7pos,X8pos]
        
        self.XAtomList = [XAtomList[i]/self.celldims for i in range(len(XAtomList))]

    #------------------------------------------------

    def initCsPositions(self) -> None:
        '''
        Guess the Cs Atomic positions.
        
        Recommend using manual init rather than this function for all subsequent relaxations. 
        
        '''
        
        X_z = self.XAtomList[4][2] # This is X5 in diagram from previous function
        Csdz = 0.9 * np.array([0,0,X_z])
        
        Cs1pos = np.array([0.5,0,self.Z]) * self.celldims + Csdz
        Cs2pos = np.array([0.5,0,self.Z]) * self.celldims - Csdz
        Cs3pos = np.array([0,0.5,self.Z]) * self.celldims + Csdz
        Cs4pos = np.array([0,0.5,self.Z]) * self.celldims - Csdz
        
        self.CsAtomList = [Cs/self.celldims for Cs in [Cs1pos,Cs2pos,Cs3pos,Cs4pos]]
        


    def readCsPositions(self,inputCsStr: str) -> None:
        '''
        Reads an input array of Cs atoms. 
        '''

        self.CsAtomList = []
        
        for line in inputCsStr.strip().split('\n'):
            elements = line.split()

            if elements:
        
                CsPos = np.array([float(elements[1]),float(elements[2]),float(elements[3])])
            
                self.CsAtomList.append(CsPos)
            
            
    #------------------------------------------------
    

    def writeCIF(self) -> None:
        '''
        
        Writes a CIF file to inspect the result in VESTA. 
        
        '''
        
    
        data = f'''data_structure
_symmetry_space_group_name_H-M    'P 1'
_cell_length_a                    {self.celldims[0]}
_cell_length_b                    {self.celldims[1]}
_cell_length_c                    {self.celldims[2]} 
_cell_angle_alpha                 90
_cell_angle_beta                  90
_cell_angle_gamma                 90
_symmetry_Int_Tables_number       1

'''

        loopdata = '''loop_
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
'''

    
        with open(f"{self.outputfile}_input.cif", "w") as file:
            
            file.write(data)
            file.write(loopdata)
            
            for i in range(len(self.BAtomList)):
                file.write(f'{self.B} {self.BAtomList[i][0]} {self.BAtomList[i][1]} {self.BAtomList[i][2]}\n')
            
            for i in range(len(self.XAtomList)):
                file.write(f'{self.X} {self.XAtomList[i][0]} {self.XAtomList[i][1]} {self.XAtomList[i][2]}\n')
                
            for i in range(len(self.CsAtomList)):
                file.write(f'Cs {self.CsAtomList[i][0]} {self.CsAtomList[i][1]} {self.CsAtomList[i][2]}\n')
        
        print(f'CIF file saved to : "{self.outputfile}_input.cif"')
            
    def writeRelaxationInput(self, FullRel: bool, nbnd: float=120) -> None:     
        
        control = f'''&control
    calculation = 'vc-relax',
    prefix = '{self.B}{self.X}_{self.beta*180/np.pi:.1f}_{self.delta*180/np.pi:.1f}',
    outdir = '{self.espressoOutdir}',
    pseudo_dir = './pseudo/',
/
    '''
    
        if FullRel == True:
            isrel = 'true'
            RelTag = '_FR'
        else:
            isrel = 'false'
            RelTag = ''
    
        system = f'''&system
    ibrav = 0,
    nat = 14,
    ntyp = 3,
    nbnd = {nbnd}
    ecutwfc = 80,
    noncolin = .{isrel}.,
    lspinorb = .{isrel}.,
    ecutrho = 320,
    occupations = 'fixed',
/
    '''
        
        electrons_ions = '''&electrons
/

&ions

/
'''

        cell = '''&cell
    cell_dofree = 'xy'
/
'''

        atomicSpecies = f'''ATOMIC_SPECIES
  {self.B} {self.Bmass} {self.B}{RelTag}.upf   
  {self.X} {self.Xmass} {self.X}{RelTag}.upf 
  {self.A} {self.Amass} {self.A}{RelTag}.upf
'''
        
    

        anstronglabel = r'{angstrom}'
        cellParameters = f'''CELL_PARAMETERS {anstronglabel}
 {self.celldims[0]} 0.000000 0.000000
 0.000000 {self.celldims[1]} 0.000000
 0.000000 0.000000 {self.celldims[2]}
        '''
        
        
        Kpoints = '''K_POINTS automatic
  5 5 1 0 0 0
    '''
        
        
        with open(f"{self.outputfile}_relax.in", "w") as file:
            
            file.write(control)
            file.write('\n')
            file.write(system)
            file.write('\n')
            file.write(electrons_ions)
            file.write('\n')
            file.write(cell)
            file.write('\n')
            file.write(atomicSpecies)
            file.write('\n')
            file.write(cellParameters)
            file.write('\n')
            
            
            file.write('ATOMIC_POSITIONS {crystal} \n')
            for i in range(len(self.BAtomList)):
                file.write(f'{self.B} {self.BAtomList[i][0]} {self.BAtomList[i][1]} {self.BAtomList[i][2]} 0 0 0 \n')
            
            for i in range(len(self.XAtomList)):
                file.write(f'{self.X} {self.XAtomList[i][0]} {self.XAtomList[i][1]} {self.XAtomList[i][2]} 0 0 0 \n')
                
            for i in range(len(self.CsAtomList)):
                file.write(f'Cs {self.CsAtomList[i][0]} {self.CsAtomList[i][1]} {self.CsAtomList[i][2]} 1 1 1 \n')
            
            file.write('\n')
            file.write(Kpoints)
            
            
        print(f'Relaxation input saved to : "{self.outputfile}_relax.in"')
        #------------------------------------------------    

'''

# Standard run (Example)
#------------------------------------------------

#Input params (Example)

celldims = np.array([  9.515504294  , 8.289755694, 22])
beta = 20
delta = 20

FullRel = True

CsAtomStr = ''
Cs               0.6072076584       -0.0456721997         0.187
Cs               0.1072076584        0.5456721997        0.187
Cs               0.3927923416        0.0456721997        -0.187
Cs              -0.1072076584        0.4543278003        -0.187
''


nbnd = 140
 
 
#cell = BuildUnitCell(celldims,beta,delta)
#cell.locateBXAtoms()
#cell.readCsPositions(CsAtomStr)
#cell.writeCIF()
#cell.writeRelaxationInput(FullRel,nbnd)

#------------------------------------------------
 
'''

