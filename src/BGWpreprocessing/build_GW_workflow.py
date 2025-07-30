'''
Builds a GW workflow for a given input 2D structure. 

Used for our study as we want to study multiple structures with multiple distortion angles

This workflow is modeleled from https://github.com/BerkeleyGW/BerkeleyGW-examples/tree/master/DFT/MoS2_SOC

William Davie 07/24
'''

# Example Atomic position list:
import numpy as np
from typing import List
import re
import os,sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))
from pathlib import Path



# these can be called using QEpreprocessing.write_inputs 
celldims = np.array([6.29669,6.29669,22.0])
AtomicPosExample='''Pb        0.0000000000        0.0000000000        0.5000000000
I         0.0000000000        0.5000000000        0.5000000000
I         0.5000000000        0.0000000000        0.5000000000
I         0.0000000000        0.0000000000        0.6442813300
I         0.0000000000        0.0000000000        0.3557186700
Cs        0.5000000000        0.5000000000        0.6200669517
Cs        0.5000000000        0.5000000000        0.3799330483'''
NoccupiedBandsExample = 60 

class buildGWworkflow():
    
    def __init__(self, celldims: np.ndarray, atomicPositions: str, NoccupiedBands: int, outdir: str, prefix: str,
                 elements: List[str]=['Pb','I','Cs'], atomicMasses: List[float]=[207.2,126.9,132.9],baseNkpoints: int=6):
        
        self.celldims = celldims
        self.atomicPositions = atomicPositions
        self.NoccupiedBands = NoccupiedBands
        
        self.elements = elements
        self.atomicMasses = atomicMasses
        self.outdir = outdir
        Path(outdir).mkdir(exist_ok=True)
        
        self.baseNkpoints = baseNkpoints
        
        self.prefix = prefix
        
        self.numAtoms = len(atomicPositions.strip().splitlines())


    def generateKgridInput(self, outputstr: str, Nkpoints: int, shiftline: str='0.0 0.0 0.0\n'):
        
        '''
        Input: Atomic positions as output by writeQEinput
        
        Output: All kgrid.x input files for:
        
        WFN WFNq
        WFN_co WFNq_co
        WFN_fi WFNq_fi
        
        self.baseNkpoints is the number of points used WFN for X X 1
        '''
        
        celldimsstr = f'''{celldims[0]:.10f} 0.0000000000 0.0000000000
0.0000000000 {celldims[1]:.10f} 0.0000000000
0.0000000000 0.0000000000 {celldims[2]:.10f}'''
        
        splitAtomicPos = self.atomicPositions.strip().splitlines()
        count = len(self.atomicPositions.strip().splitlines())
        atomstrfinal = f'{count}\n'

        for line in splitAtomicPos:
            lineSplit = re.split('[ ]+',line)
            label = self.elements.index(lineSplit[0]) + 1
            
            atomstrfinal += f'{label} {lineSplit[1]} {lineSplit[2]} {lineSplit[3]} \n'
        

        #-------------WFN-------------
        with open(outputstr,'w') as file:
        
            file.write(f'{Nkpoints} {Nkpoints} 1\n')
            file.write('0.0 0.0 0.0\n')
            file.write(shiftline)
            file.write('\n')
            file.write(celldimsstr + '\n')
            file.write(atomstrfinal)
            file.write(f'{self.baseNkpoints*5} {self.baseNkpoints*5} {self.baseNkpoints*5*5}\n')
            file.write('.false.')
            
    
    def kgridWorkflow(self, basic: bool=True):
        '''
        Calls generate k point grid for all 
        '''
        
        kgridDir = self.outdir+'/00-kgrids/'
        Path(kgridDir).mkdir(exist_ok=True)
        
        
        # For a 2D structure we using subsampling to optain WFNq, hence why not included here
        
        #-------------WFN------------- 
        if basic == False:
            self.generateKgridInput(kgridDir+'WFN_kgrid.inp',self.baseNkpoints)
        #-------------WFN_co------------- 
        self.generateKgridInput(kgridDir+'/WFN_co_kgrid.inp',self.baseNkpoints*2)
        #-------------WFNq_co------------- 
        self.generateKgridInput(kgridDir+'WFNq_co_kgrid.inp',self.baseNkpoints*2,shiftline='0.001 0.0 0.0\n')
        #-------------WFN_fi------------- 
        self.generateKgridInput(kgridDir+'WFN_fi_kgrid.inp',self.baseNkpoints*6)
        #-------------WFNq_fi------------- 
        self.generateKgridInput(kgridDir+'WFNq_fi_kgrid.inp',self.baseNkpoints*6,shiftline=f'0.001 0.0 0.0\n')
        
    def change2uniformKgrid(self,kgridOutFile: str, kgridLogFile: str):
        '''
        use when want to prevent BZ folding for the absorption calculation. 

        only really used for WFN_fi
        
        allows us to use use_velocity.
        '''
        
        logfile = open(self.outdir + '/00-kgrids/' + kgridLogFile,'r')
        lines = logfile.readlines()
        
        for i, line in enumerate(lines):
            
            if 'k-points in the original uniform grid' in line:
                
                startIndex = i 
                numKpoints = int(lines[i+1])
                
        output = f'''K_POINTS crystal
    {numKpoints}\n'''
        
        for i, line in enumerate(lines[startIndex+2:startIndex+2+numKpoints]):
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
            output += f'{linesplit[1]} {linesplit[2]} {linesplit[3]} 1.0\n'
            
        with open(self.outdir + '/00-kgrids/' +kgridOutFile,'w') as file:
            file.write(output)
            
        return None

    
    def writeNSCFinput(self, nbnd: int, postfix: str, kgridfile: str='',calculation: str='nscf', outdir: str=''):
        
        Path(outdir).mkdir(exist_ok=True)
        
        '''
        
        '''
        # <!> some repition from writeQEinput could be avoided with smarter workflow but for now we re-write.
        
        if calculation=='scf': wfcolStr = '.true.'
        else: wfcolStr = '.false.'
    
        
        control = f'''&control
    calculation = '{calculation}',
    prefix = '{self.prefix}',
    outdir = './',
    disk_io = 'low',
    pseudo_dir = '../pseudo',
    verbosity = 'low',
    wf_collect = {wfcolStr}
/'''
        system = f'''&system
    ibrav = 0,
    nat = {self.numAtoms},
    ntyp = 3,
    nbnd = {nbnd},
    ecutwfc = 90,
    ecutrho = 360,
    noncolin = .true.,
    lspinorb = .true.,
    occupations = 'fixed',
    assume_isolated = '2D'
/'''
        electrons = '''&electrons
    conv_thr = 1.0d-8,
/'''
    
    # could be edited for more elements with more 
        atomicSpecies = f'''ATOMIC_SPECIES 
  {self.elements[0]} {self.atomicMasses[0]} {self.elements[0]}_FR.upf   
  {self.elements[1]} {self.atomicMasses[1]} {self.elements[1]}_FR.upf 
  {self.elements[2]} {self.atomicMasses[2]} {self.elements[2]}_FR.upf
'''

        CellParams = f'''CELL_PARAMETERS (angstrom)
   {celldims[0]}   0.000000000   0.000000000
   0.000000000   {celldims[1]}   0.000000000
   0.000000000   0.000000000   {celldims[2]}
'''
            
        AtomicPos = f'''ATOMIC_POSITIONS (crystal)
{self.atomicPositions}
            '''
            
        if calculation == 'scf':
            Kpoints = f'''K_POINTS automatic
{self.baseNkpoints*2} {self.baseNkpoints*2} 1 0 0 0
'''
        else:
            try:
                kgridlines = open(kgridfile,'r').readlines()
                kgrid = ''.join(kgridlines)
            
                Kpoints = f'''{kgrid}
                '''
            except: # this is specifically for the case WFNq where the k points are initially unknown. 
                print('<!> Warning K point file not found, check if input is for 02-wfn.')
                Kpoints = f'''K_POINTS crystal 
<!> Not set
'''
                
    
        with open(f"{outdir}/{self.prefix}_{postfix}.in", "w") as file:
        
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
            
    def writePW2BGWinput(self, Nkpoints: int, wfnFilename: str='WFN',
                         rhoFile: bool=True, vxcFile: bool=True, shiftarray: np.ndarray=np.array([0,0,0]), 
                         outdir: str=''):
        
        if rhoFile: rhoFileStr = '.true.'
        else: rhoFileStr = '.false.'

        if vxcFile: vxcFileStr = '.true.'
        else: vxcFileStr = '.false.'

        # WARNING <!> check diag and off_diag elements.   
            
        input=f'''&input_pw2bgw
  prefix = '{self.prefix}'
  real_or_complex = 2 
  wfng_flag = .true. 
  wfng_file = '{wfnFilename}'
  wfng_kgrid = .true.
  wfng_nk1 = {Nkpoints}
  wfng_nk2 = {Nkpoints}
  wfng_nk3 = 1
  wfng_dk1 = {shiftarray[0]}
  wfng_dk2 = {shiftarray[1]}
  wfng_dk3 = {shiftarray[2]}
  wfng_occupation = .false.
  wfng_nvmin = 0
  wfng_nvmax = 0
  rhog_flag = {rhoFileStr}
  rhog_file = 'RHO'
  vxc_flag = {vxcFileStr}
  vxc_file = 'vxc.dat'
  vxc_diag_nmin = 1
  vxc_diag_nmax = 12
 /'''
 
        with open(f"{outdir}/pw2bgw.in", "w") as file:
            
            file.write(input)
            
        return None
           

            
    def buildQEworkflow(self, NunoccupiedBands: int, basic: bool=True):
        '''
        Generates all pw.x input files calculations to 
        '''
        
        kgridDir = self.outdir + '/00-kgrids/'
        
        #----01-SCF------
        
        self.writeNSCFinput(self.NoccupiedBands+10,'scf',calculation='scf',outdir=self.outdir+'/01-scf')
        self.writePW2BGWinput(self.baseNkpoints*2,wfnFilename='WFN',rhoFile=True,vxcFile=False,outdir=self.outdir+'/01-scf')
        
        if basic == False: # basic = shortest workflow to obtain exciton energies. 
         
            #----02-WFN------
            self.writeNSCFinput(self.NoccupiedBands+NunoccupiedBands,'wfn',kgridfile=kgridDir+'WFNq_kgrid.out') # large number of unoccupied bands for convergence. 
            
            #----03-WFNq-----
            #self.writeNSCFinput(self.NoccupiedBands+10,'wfnq')
            
            #----08-bands----
            #self.writeNSCFinput(self.NoccupiedBands+NunoccupiedBands)
             
        
        #----04-WFN_co---  #WARNING <!> need to change diagonalization here. 
        
        self.writeNSCFinput(self.NoccupiedBands+NunoccupiedBands,'wfn_co',kgridfile=kgridDir+'WFN_co_kgrid.out',outdir=self.outdir+'/04-wfn_co')  # 
        self.writePW2BGWinput(self.baseNkpoints*2,wfnFilename='WFN',rhoFile=False,vxcFile=True,outdir=self.outdir+'/04-wfn_co')
        
        #----05-WFNq_co--
        
        self.writeNSCFinput(self.NoccupiedBands+10,'wfnq_co',kgridfile=kgridDir+'WFNq_co_kgrid.out',outdir=self.outdir+'/05-wfnq_co') 
        self.writePW2BGWinput(self.baseNkpoints*2,wfnFilename='WFNq',rhoFile=False,vxcFile=False,outdir=self.outdir+'/05-wfnq_co',shiftarray=np.array([0.001,0,0]))
        
        #----06-WFN_fi---
        self.change2uniformKgrid('WFN_fi_kgrid.out',"WFN_fi_kgrid.log")
        self.writeNSCFinput(self.NoccupiedBands+10,'wfn_fi',kgridfile=kgridDir+'WFN_fi_kgrid.out',outdir=self.outdir+'/06-wfn_fi') # make sure to set this to the unfolded BZ
        self.writePW2BGWinput(self.baseNkpoints*6,wfnFilename='WFN_fi',rhoFile=False,vxcFile=False,outdir=self.outdir+'/06-wfn_fi',shiftarray=np.array([0,0,0])) # 1/(2*self.baseNkpoints*6)
        
        #----07-WFNq_fi--
        
        self.writeNSCFinput(self.NoccupiedBands+10,'wfnq_fi',kgridfile=kgridDir+'WFNq_fi_kgrid.out',outdir=self.outdir+'/07-wfnq_fi') 
        self.writePW2BGWinput(self.baseNkpoints*6,wfnFilename='WFNq_fi',rhoFile=False,vxcFile=False,outdir=self.outdir+'/07-wfnq_fi',shiftarray=np.array([0.001,0,0]))
        
        
    def k2epsilonQ(self, inputKpath: str, shift: float=0.001) -> str:
        '''
        Given a set of k points return a qpoint grid for epsilon:
        '''

        kpointfile = open(inputKpath,'r')
        lines = kpointfile.readlines()
        
        outputQpoints = ''
        
        for line in lines[2:]:
            
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
            kpoint = np.array([linesplit[0],linesplit[1],linesplit[2]],dtype=np.float32)
            
            # add shift to prevent 1/|q| = inf
            if kpoint[0] == kpoint[1] == kpoint[2] == 0:
                kpoint[0] += shift
                endstr = '1.0 1'
            else:
                endstr = '1.0 0'
                
            outputQpoints += f'{kpoint[0]:.10f} {kpoint[1]:.10f} {kpoint[2]:.10f} {endstr}\n'
                
        return outputQpoints
        
        

    def writeEpsilonIn(self, qpoints: str, outdir: str):
        
        Path(outdir).mkdir(exist_ok=True)
        
        cutoff = 5.0
        nbnds = 100
        
        input = f'''# epsilon.inp
epsilon_cutoff           {cutoff}
number_bands             {nbnds}
cell_slab_truncation
frequency_dependence 0
degeneracy_check_override
begin qpoints
{qpoints}end
        '''
        
        with open(outdir+'/epsilon.inp','w') as file:
            
            file.write(input)
            
        return None 
    
    def writeSigmaIn(self, inputKpath: str, outdir: str, subsample: bool=False):
        
        Path(outdir).mkdir(exist_ok=True)
        
        # use inputKpath to input k points - could be optimised with k2epsilonQ but leaving for now 
        
        kpoints = ''
        
        klines = open(inputKpath,'r').readlines()
        
        for line in klines[2:]:
            
            linestrip = line.strip()
            linesplit = re.split('[ ]+',linestrip)
            
            kpoint = np.array([linesplit[0],linesplit[1],linesplit[2]],dtype=np.float32)
            
            kpoints += f'{kpoint[0]:.10f} {kpoint[1]:.10f} {kpoint[2]:.10f} 1.0\n'
        
        
        if subsample==True: subsampleStr = '\nsubsample\n'
        else: subsampleStr = ''
        
        # WARNING <!> these band indicies do not work, must be paired with degenerate parters. 
        input = f'''# sigma.inp
band_index_min   {self.NoccupiedBands-12} 
band_index_max   {self.NoccupiedBands+12}
screening_semiconductor
degeneracy_check_override
cell_slab_truncation
exact_static_ch 1{subsampleStr}
begin kpoints
{kpoints}end'''
            
        with open(outdir+'/sigma.inp','w') as file:
            file.write(input)
            
        return None
    
    
    def writeKernelIn(self, outdir: str):
        
        Path(outdir).mkdir(exist_ok=True)
        
        nVbnds = 8 # need converging 
        nCbnds = 8 
        
        input = f'''number_val_bands {nVbnds}
number_cond_bands {nCbnds}
cell_slab_truncation
use_symmetries_coarse_grid
low_comm
screening_semiconductor'''

        with open(outdir + '/kernel.inp','w') as file:
            file.write(input)
            
        return None 
    
    
    def writeAbsorptionIn(self, outdir: str):
        
        Path(outdir).mkdir(exist_ok=True)
        
        nVband_co = 8 # same as kernel
        nCband_co = 8 
        
        nVband_fi = 2
        nCband_fi = 2
        
        input = f'''spinor
diagonalization

number_val_bands_coarse {nVband_co}
number_val_bands_fine {nVband_fi}
number_cond_bands_coarse {nCband_co}
number_cond_bands_fine {nCband_fi}

use_symmetries_coarse_grid
use_symmetries_fine_grid
use_symmetries_shifted_grid

screening_semiconductor
cell_slab_truncation
use_velocity
gaussian_broadening
energy_resolution 0.02
eqp_co_corrections
        '''
        
        with open(outdir + '/absorption.inp','w') as file:
            
            file.write(input)
            
        return None
            
        
        
        
    def buildGWworkflow(self):
        
        
        #----EpsilonCo-----
        epsilon_co_qpoints = self.k2epsilonQ(self.outdir+'/00-kgrids/WFN_co_kgrid.out')
        self.writeEpsilonIn(qpoints=epsilon_co_qpoints,outdir=self.outdir+'/11-epsilon_co')
        
        #-------Sigma------
        self.writeSigmaIn(self.outdir+'/00-kgrids/WFN_co_kgrid.out',self.outdir+'/12-sigma')
        
    
        #------Kernel------
        self.writeKernelIn(self.outdir+'/16-kernel') 
        
        #----Absorption----
        
        self.writeAbsorptionIn(self.outdir+'/17-absorption')
        
