'''
A script to plot the projected density of states for an undistored CsPbI structure:

used in my report to highlight where which states contribute to the band structure.

William Davie 07/25
'''

import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..")))
from QEpostprocessing.read_QE_output import readQEouput,sumDirDOS
from plotting_functions import runFormat,fetchColour

# directory

runFormat()


DOSdir = './src/scripts/scriptdata/projwfc_0.0_0.0'

#-------------------Total DOS-----------------------

dataTotalDOS = readQEouput(DOSdir + '/' + 'CsPbI_0.0_0.0_proj.pdos_tot')
Evals,totalDOS = dataTotalDOS.readProjDOS('CsPbI_0.0_0.0_proj.pdos_tot')

Evals,CsDOS = sumDirDOS(directory=DOSdir,AtomTag='Cs',OrbitalTag='')
Evals,Pb_p_DOS = sumDirDOS(directory=DOSdir,AtomTag='Pb',OrbitalTag='p')
Evals,Pb_s_DOS = sumDirDOS(directory=DOSdir,AtomTag='Pb',OrbitalTag='s')
Evals,I_p_DOS = sumDirDOS(directory=DOSdir,AtomTag='I',OrbitalTag='p')

Evals = Evals + 4.1688
start=0


fig,ax = plt.subplots()
ax.set(xlabel='E (eV)')
ax.plot(Evals[start:],totalDOS[start:],color='black',label='Total')
ax.margins(x=0)
ax.set_yticks([])
# block1 Cs - I - Pb
block1end = 500
ax.fill_between(Evals[start:block1end],CsDOS[start:block1end],label='$\mathrm{Cs}^{+}$',color='grey',alpha=0.8)
ax.fill_between(Evals[start:block1end],I_p_DOS[start:block1end],label='$\mathrm{I}^{-}$: $p$',color=fetchColour('lightred'),alpha=0.8)
ax.fill_between(Evals[start:block1end],Pb_p_DOS[start:block1end],label='$\mathrm{Pb}^{2+}$: $p$',color=fetchColour('blue'),alpha=0.8)
# block2 I - Pb - CS
block2end = 1100
ax.fill_between(Evals[block1end:block2end],I_p_DOS[block1end:block2end],color=fetchColour('lightred'),alpha=0.8)
ax.fill_between(Evals[block1end:block2end],Pb_p_DOS[block1end:block2end],color=fetchColour('blue'),alpha=0.8)
ax.fill_between(Evals[block1end:block2end],CsDOS[block1end:block2end],color='grey',alpha=0.8)
# block3 Cs - Pb - I 
ax.fill_between(Evals[block2end:],CsDOS[block2end:],color='grey',alpha=0.8)
ax.fill_between(Evals[block2end:],Pb_p_DOS[block2end:],color=fetchColour('blue'),alpha=0.8)
ax.fill_between(Evals[block2end:],I_p_DOS[block2end:],color=fetchColour('lightred'),alpha=0.8)
ax.set(xlim=[-10,5])
ax.legend()
plt.savefig('projDOS.pdf')
plt.show()

