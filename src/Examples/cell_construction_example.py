'''
Example of building a 2D pervoskite unit cell
'''

from QEpreprocessing.build_unit_cell import BuildUnitCell
import numpy as np

celldims = np.array([    9.651083020 , 8.032582030, 22])
beta = 0
delta = 0

FullRel = True

CsAtomStr = '''
Cs               0.6072076584       -0.0456721997        0.1872100505
Cs               0.1072076584        0.5456721997        0.1872100505
Cs               0.3927923416        0.0456721997       -0.1872100505
Cs              -0.1072076584        0.4543278003       -0.1872100505
'''

nbnd = 140

cell = BuildUnitCell(celldims,beta,delta)
cell.locateBXAtoms()
cell.readCsPositions(CsAtomStr)
cell.writeCIF()
cell.writeRelaxationInput(FullRel,nbnd)