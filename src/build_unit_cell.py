

'''
Constructs a `double model' unit cell for a distorted 2D Hallide pervoskite with a caesium cation.

General formula: (Cs)_2 (B X_4)

Note that this cell is build via emperical approximations and should be relaxed via DFT. 

William Davie 17/04/25
'''

import numpy as np
from typing import List

#Input params (Example)

B = 'Pb' 
X = 'I'  
celldims = np.array([8.977963,8.977963,16])
beta = 0
delta = 0 


CsAtomStr = '''
Cs               0.5000000000        0.0000000000        0.1608654327
Cs               0.0000000000        0.5000000000        0.1608564732
Cs               0.5000000000        0.0000000000       -0.1608654327
Cs               0.0000000000        0.5000000000       -0.1608564732
'''

AtomicSpecies = ''' 
ATOMIC_SPECIES
  Pb 207.2 Pb.pbe-dn-kjpaw_psl.1.0.0.UPF
  I  126.9 I.pbe-n-kjpaw_psl.1.0.0.UPF
  Cs 132.9 Cs.pbe-spn-kjpaw_psl.1.0.0.UPF
'''



class BuildUnitCell():
    
    def __init__(self, B: str, X: str, celldims: np.ndarray, beta: float, delta: float):
        
        self.B = B
        self.X = X
        
        self.beta = beta*np.pi/180
        self.delta = delta*np.pi/180
        
        self.celldims = celldims
        
        self.filename = f'Cs{B}{X}_{beta}_{delta}'
        
        self.outdir = '/local/data/public/wd324/QEout/' # Change when needed
        
        self.XAtomList = []
        self.BAtomList = []
        self.CsAtomList = []
        

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
        
        B1pos = np.array([0.0,0.0,0.0])*self.celldims
        B2pos = np.array([0.5,0.5,0.0])*self.celldims
        
        # in plane X atom calculation:
        
        # initial positions:
        
        X1pos = np.array([0.25,0.25,0.0])*self.celldims
        X2pos = np.array([0.75,0.75,0.0])*self.celldims
        X3pos = np.array([0.25,0.75,0.0])*self.celldims
        X4pos = np.array([0.75,0.25,0.0])*self.celldims

        # B - B bond length 

        B_BLength = np.sqrt(self.celldims[0]**2 + self.celldims[1]**2)
        
        cosTheta = (self.celldims[1]/2) / (B_BLength/2) # could also use celldims[0] but easier option is to just switch 0 and 1 within input


        # based on the constraints of delta and beta on the Pb - I - Pb bond, we arrive at:
        
        L = np.sqrt( ( (1/4 * B_BLength**2)*(1 + np.tan(self.beta)**2) ) / ( cosTheta**2 * np.tan(self.delta)**2 + 1 )) / 2 

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
        
        Cs1pos = np.array([0.5,0,0]) * self.celldims + Csdz
        Cs2pos = np.array([0.5,0,0]) * self.celldims - Csdz
        Cs3pos = np.array([0,0.5,0]) * self.celldims + Csdz
        Cs4pos = np.array([0,0.5,0]) * self.celldims - Csdz
        
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

    
        with open(f"{self.filename}_input.cif", "w") as file:
            
            file.write(data)
            file.write(loopdata)
            
            for i in range(len(self.BAtomList)):
                file.write(f'{self.B} {self.BAtomList[i][0]} {self.BAtomList[i][1]} {self.BAtomList[i][2]}\n')
            
            for i in range(len(self.XAtomList)):
                file.write(f'{self.X} {self.XAtomList[i][0]} {self.XAtomList[i][1]} {self.XAtomList[i][2]}\n')
                
            for i in range(len(self.CsAtomList)):
                file.write(f'Cs {self.CsAtomList[i][0]} {self.CsAtomList[i][1]} {self.CsAtomList[i][2]}\n')
            
            
    def writeRelaxationInput(self, AtomicSpeciesStr: str) -> None:     
        
        control = f'''&control
    calculation = 'vc-relax',
    prefix = '{self.B}{self.X}_{self.beta}_{self.delta}',
    outdir = '/local/data/public/wd324/QEout/',
    pseudo_dir = './pseudo/',
    forc_conv_thr = 1.0d-4,
/
    '''

        system = '''&system
    ibrav = 0,
    nat = 14,
    ntyp = 3,
    ecutwfc = 80,
    ecutrho = 320,
    occupations = 'fixed',
/
    '''
        
        electrons_ions = '''&electrons
    conv_thr = 1.0d-8,
/

&ions

/
'''

        cell = '''&cell
    cell_dofree = 'xy'
/
'''

        atomic_species = AtomicSpeciesStr
        
    

        anstronglabel = r'{angstrom}'
        cell_parameters = f'''CELL_PARAMETERS {anstronglabel}
 {self.celldims[0]} 0.000000 0.000000
 0.000000 {self.celldims[1]} 0.000000
 0.000000 0.000000 {self.celldims[2]}
        '''
        
        
        Kpoints = '''K_POINTS automatic
  5 5 1 0 0 0
    '''
        
        
        with open(f"{self.filename}_relax.in", "w") as file:
            
            file.write(control)
            file.write('\n')
            file.write(system)
            file.write('\n')
            file.write(electrons_ions)
            file.write('\n')
            file.write(cell)
            file.write('\n')
            file.write(atomic_species)
            file.write('\n')
            file.write(cell_parameters)
            file.write('\n')
            
            
            file.write('ATOMIC_POSITIONS {crystal} \n')
            for i in range(len(self.BAtomList)):
                file.write(f'{self.B} {self.BAtomList[i][0]} {self.BAtomList[i][1]} {self.BAtomList[i][2]}\n')
            
            for i in range(len(self.XAtomList)):
                file.write(f'{self.X} {self.XAtomList[i][0]} {self.XAtomList[i][1]} {self.XAtomList[i][2]}\n')
                
            for i in range(len(self.CsAtomList)):
                file.write(f'Cs {self.CsAtomList[i][0]} {self.CsAtomList[i][1]} {self.CsAtomList[i][2]}\n')
            
            file.write('\n')
            file.write(Kpoints)
            
            
        #------------------------------------------------
        
        


# Standard run. 
#------------------------------------------------
 
 
cell = BuildUnitCell(B,X,celldims,beta,delta)
cell.locateBXAtoms()
cell.readCsPositions(CsAtomStr)
cell.writeCIF()
cell.writeRelaxationInput(AtomicSpecies)
