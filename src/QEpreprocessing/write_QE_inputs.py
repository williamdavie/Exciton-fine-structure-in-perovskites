'''

Given a CIF file of structures, write QE inputs (SCF, NSCF, etc).

General formula: (Cs)_2 (B X_4)

Assumes the following formatting title as: data_beta_delta

'''

import numpy as np
from typing import List
import math, re, sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))
from perovskite_2D import perovskite2D
from QEpostprocessing.read_QE_output import readQEouput

#------------------------------------------------

class writeQEinput(perovskite2D):
    
    def __init__(self, filename: str, FullRel=True, 
                 startline: int=0, endline: int=None) -> None:
        
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
        self.Kpath = f'K_POINTS crystal\n{(numSteps*(numHSpoints-1)+1)}\n'
        
        for i in range(1,numHSpoints):
            
            start = highSymPoints[i-1]
            end = highSymPoints[i]
            
            for j in range(numSteps):
                
                point = start + (end-start) * j/numSteps
                
                self.Kpath += f'{point[0]:.10f} {point[1]:.10f} {point[2]:.10f} 1\n'
        
        # add the last point
        finalpoint = highSymPoints[len(highSymPoints)-1]
        self.Kpath += f'{finalpoint[0]:.10f} {finalpoint[1]:.10f} {finalpoint[2]:.10f} 1\n'
        
        if outputFile == True:
            with open(f'kpoints_{highSymPoints[0]}-{highSymPoints[numHSpoints]}.in', 'w') as f:
                f.write(self.Kpath)
    
    #---------------------------------------------
    
    def defineKmesh(self, Nkx: int, Nky: int, Nkz: int, scale: float=1/3 ,outputFile: bool=False):
        '''
        Defines a gamma centred k point mesh, used for effective mass calculations
        '''
        assert Nkx % 2 == 1 and Nky % 2 == 1 and Nkz % 2 == 1, 'Must be odd to be gamma centred'
        
        self.Kpath = f'K_POINTS crystal\n{Nkx*Nky*Nkz}\n'
        Karray = np.zeros((Nkx,Nky,Nkz,3))
        for i in range(Nkx):
            for j in range(Nky):
                for k in range(Nkz):
                    
                    # -math.floor(N/2)/N centres the grid around Gamma as long as the first assert is satisfied.
                    
                    x = (i / Nkx - math.floor(Nkx/2)/Nkx) * scale
                    y = (j / Nky - math.floor(Nky/2)/Nkx) * scale
                    z = (k / Nkz - math.floor(Nkz/2)/Nkx) * scale
                    
                    self.Kpath += f'{x:.10f} {y:.10f} {z:.10f} 1\n'
                    Karray[i,j,k] = np.array([x,y,z])
                

        if outputFile == True:
            with open(f'kpoints_{Nkx}x{Nky}x{Nkz}_gamma.in', 'w') as f:
                f.write(self.Kpath)
            np.save(f'kpoints_{Nkx}x{Nky}x{Nkz}_gamma',Karray.reshape(Nkx*Nky*Nkz,3))
                
            
        
    #---------------------------------------------
        

    def writePWSCF(self, calculation: str, outdir: str, 
                   nbnd: int=140, assume_isolated: str='2D', verbosity: str='low',
                   postfix: str='') -> None:
        '''
        Writes a PWSCF input for each structure.
        '''
    
        for i in range(self.numStructures):
            
            control = f'''&control
    calculation = '{calculation}',
    prefix = '{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}',
    outdir = '/home/wd324/rds/hpc-work/CsPbI_Out',
    pseudo_dir = '/home/wd324/pseudo',
    verbosity = '{verbosity}'
/'''
            system = f'''&system
    ibrav = 0,
    nat = 14,
    ntyp = 3,
    nbnd = {nbnd},
    ecutwfc = 90,
    ecutrho = 360,
    noncolin = .{self.isRel}.,
    lspinorb = .{self.isRel}.,
    occupations = 'fixed',
    assume_isolated = '{assume_isolated}'
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
    
            with open(f"./{outdir}/{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}_{calculation}{postfix}.in", "w") as file:
        
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
    prefix = '{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}',
    outdir = '/local/data/public/wd324/QEout/',
    filband='{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}_band.dat'
/  
    '''
            with open(f"./{outdir}/{self.prefix}_{self.beta_vals[i]}_{self.delta_vals[i]}_bands.in", "w") as file:
         
                file.write(input)