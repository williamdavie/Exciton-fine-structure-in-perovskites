'''

Given a CIF file of structures, write QE inputs (SCF, NSCF, etc).

General formula: (Cs)_2 (B X_4)

Assumes the following formatting title as: data_beta_delta

'''

import numpy as np
import re
from typing import List
from perovskite_2D import perovskite2D

#------------------------------------------------
# Input 

filename = '(Cs)2_PbI4_structures_0_2.5.cif'

Kpath = [
    np.array([0.5000, 0.5000, 0.0]),
    np.array([0.0000, 0.0000, 0.0]),
    np.array([0.5000, 0.0000, 0.0])
]

#------------------------------------------------

class buildInputFiles(perovskite2D):
    
    def __init__(self, filename: str, FullRel=True, startline: int=0, endline: int=None) -> None:
        
        super().__init__(B='Pb', Bmass=207.2, X='I', Xmass='126.9')
        
        self.prefix = f'{self.A}{self.B}{self.X}'

        file = open(filename,'r')
        self.lines = file.readlines()
        
        # data to collect
        
        self.beta_vals = []
        self.delta_vals = []
        
        self.a_vals = []
        self.b_vals = []
        self.c_vals = []
        self.ATOMIC_POSITIONS_list = []
        
        # check if endline specified
        if endline == None:
            endline = len(self.lines)
        
        for i,line in enumerate(self.lines[startline:endline]):
            
            line.strip()
            lineSplit = re.split('[ ]+',line)
             
             # beta and delta
            if '#' in line:
                continue
             
            if 'data_' in line:
                print('Line ',line)
                underscore = [i for i, underscore in enumerate(line) if underscore == '_']
                self.beta_vals.append(float(line[underscore[0]+1:underscore[1]]))
                self.delta_vals.append(float(line[underscore[1]+1:len(line)+1]))
                

            # cell dims
            
            if '_cell_length_a' in line:
                self.a_vals.append(float(lineSplit[1]))
                
            if '_cell_length_b' in line:
                self.b_vals.append(float(lineSplit[1]))
                
            if '_cell_length_c' in line:
                self.c_vals.append(float(lineSplit[1]))
            
            # atomic positions 
            
            if  '_atom_site_fract_z' in line:
                
                AtomicPos = "\n".join(line.strip() for line in self.lines[i+1:i+15])
                
                self.ATOMIC_POSITIONS_list.append(AtomicPos)
                
    
        # other variables
        
        self.numStructures = len(self.beta_vals)
        self.Kpath = None

        if FullRel == True:
            self.isRel = 'true'
            self.RelTag = '_FR'
        else:
            self.isRel = 'false'
            self.RelTag = ''
            
        # used for input files, constant throughout


            
    #---------------------------------------------
    
    def defineKpath(self,
                    highSymPoints: List[np.ndarray], numSteps: int, outputFile: bool=False):
        '''
        For band structure calculations defines a path.
        
        numSteps : steps between each k point.
        
        '''
        numHSpoints = len(highSymPoints)
        self.Kpath = f'K_POINTS crystal\n{len(numSteps*(numHSpoints-1)+1)}\n'
        
        for i in range(1,numHSpoints):
            
            start = highSymPoints[i-1]
            end = highSymPoints[i]
            
            for j in range(numSteps):
                
                point = start + (end-start) * j/numSteps
                
                self.Kpath += f'{point[0]:.10f} {point[1]:.10f} {point[2]:.10f} 1\n'
        
        if outputFile == True:
            with open(f'kpoints_{highSymPoints[0]}-{highSymPoints[numHSpoints]}.in', 'w') as f:
                f.write(self.Kpath)
    
    #---------------------------------------------
    
    def defineKmesh(self, Nkx: int, Nky: int, Nkz: int, scale: float=1/3 ,outputFile: bool=False):
        '''
        Defines a gamma centred k point mesh, used for effective mass calculations
        '''
        
        self.Kpath = f'K_POINTS crystal\n{Nkx*Nky*Nkz}\n'
        for i in range(Nkx):
            for j in range(Nky):
                for k in range(Nkz):
                    
                    x = i / Nkx * scale
                    y = j / Nky * scale
                    z = k / Nkz * scale

                self.Kpath += f'{x:.10f} {y:.10f} {z:.10f} 1\n'

        if outputFile == True:
            with open(f'kpoints_{Nkx}x{Nky}x{Nkz}_gamma.in', 'w') as f:
                f.write(self.Kpath)
        
    #---------------------------------------------
        

    def writePWSCF(self, calculation: str, outdir: str) -> None:
        '''
        Writes a PWSCF input for each structure.
        '''
    
        for i in range(self.numStructures):
            
            control = f'''&control
    calculation = '{calculation}',
    prefix = '{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}',
    outdir = '/local/data/public/wd324/QEout/',
    pseudo_dir = './pseudo/',
/'''
            system = f'''&system
    ibrav = 0,
    nat = 14,
    ntyp = 3,
    nbnd = 150,
    ecutwfc = 80,
    ecutrho = 320,
    noncolin = .{self.isRel}.,
    lspinorb = .{self.isRel}.,
    occupations = 'fixed',
/'''
            electrons = '''&electrons
    conv_thr = 1.0d-8,
/'''

            atomicSpecies = f'''ATOMIC_SPECIES
  {self.B} {self.Bmass} {self.B}{self.RelTag}.upf   
  {self.X} {self.Xmass} {self.X}{self.RelTag}.upf 
  {self.A} {self.Amass} {self.A}{self.RelTag}.upf
'''

            CellParams = f'''CELL_PARAMETERS (angstrom)
   {self.a_vals[i]}   0.000000000   0.000000000
   0.000000000   {self.b_vals[i]}   0.000000000
   0.000000000   0.000000000  {self.c_vals[i]}
'''
            
            AtomicPos = f'''ATOMIC_POSITIONS (crystal)
{self.ATOMIC_POSITIONS_list[i]}
            '''
            
            if calculation == 'scf':      
                Kpoints = '''K_POINTS automatic
    5 5 1 0 0 0
        '''
            else:
                assert self.Kpath is not None, "Must define a Kpath defineKpath()"
                Kpoints = self.Kpath
    
            with open(f"./{outdir}/{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}_{calculation}.in", "w") as file:
        
                file.write(control)
                file.write('\n')
                file.write('\n')
                file.write(system)
                file.write('\n')
                file.write('\n')
                file.write(electrons)
                file.write('\n')
                file.write(atomicSpecies)
                file.write('\n')
                file.write(CellParams)
                file.write('\n')
                file.write(AtomicPos)
                file.write('\n')
                file.write(Kpoints)
                
    def writeBands(self, outdir: str) -> None:
        '''
        Writes input for bands.x
        '''
        
        for i in range(self.numStructures):
            
            input = f'''&bands
    prefix = 'CsSnI_{self.beta_vals[i]}_{self.delta_vals[i]}',
    outdir = '/local/data/public/wd324/QEout/',
    filband='CsSnI_{self.beta_vals[i]}_{self.delta_vals[i]}_band.dat'
/  
    '''
            with open(f"./{outdir}/{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}_bands.in", "w") as file:
         
                file.write(input)
            
        
        
# Standard run 
#------------------------------------------------
 
data = buildInputFiles(filename)
data.writePWSCF('scf','SCF_0_2.5_files')
data.defineKmesh(5,5,1,scale=1/2,outputFile=True)
data.writePWSCF('nscf','SCF_0_2.5_files')
        
#------------------------------------------------
 
        