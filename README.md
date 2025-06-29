
# Tuning the Exciton Fine Structure in 2D Perovskites

Documentation for the code supporting the completion of my Masters thesis.  

Department of Materials Science and Metallurgy,  
Cavendish Laboratory Department of Physics,  
University of Cambridge

### Contents

1. [Quantum Expresso Pre-Processing](#quantum-espresso-pre-processing)
2. [Figure Creation](#figure-creation)

## Quantum Espresso Pre-Processing

### pervoskite_2D.py

Defines a base class used throughout, trivial - formatting.

### build_unit_cell.py

Sets up a unit cell with general formula Cs₂BX₄ according to distortion angles β and δ:

<div align="center">
  <img src="https://github.com/williamdavie/Exciton-fine-structure-in-perovskites/blob/main/figures/cell_marked_angles.png" width="33%">
</div>

Based on emperical investigation, the cell is setup to satisfy the following constraints:

- Angles β and δ are fixed to a pre-set value.
- The bond length L (B-X) is constant throughout the cell.
- The angle between the in-plane and out-of-plane X-site atoms is 90°:
<div align="center">
  <img src="https://github.com/williamdavie/Exciton-fine-structure-in-perovskites/blob/main/figures/90degree_marked.png" width="20%">
</div>

Employing these constraints leads to the following restriction on the bond length:

$$L = \sqrt{\frac{L_0^2(1 + \tan{\beta}^2)}{(1 + \cos{\theta}^2\tan{\delta}^2)}}$$

where $L_0 = \frac{1}{4}\sqrt{a^2 + b^2}$ (see report for detail).

Once the cell is constructed it may be viewed via an output CIF file **writeCIF(self)** and relaxed via QE **writeRelaxationInput(self)**.

Relaxation details specific to our work:

- 'vc-relax' with largely default input params [QE](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
- Force convergence threshold of 1e-3 (a.u.).
- Norm conserving fully relativistic pseudo potentials NC FR [pseudo-dojo](https://www.pseudo-dojo.org/) NC FR (ONCVPSP v0.4).
- Spin orbit coupling included.

All final structures in CIF format are given in [./Structures](https://github.com/williamdavie/Exciton-fine-structure-in-perovskites/tree/main/structures).
Structures are considered `final' when force convergence (1e-3) is reached with angles within 2.d.p of the desired value. 

### write_inputs.py

Given a set of structures in a CIF file quantum espresso input files can be constructed, including .in for:

- SCF
- NSCF, K points for evaluation are calculated as:
    - A gamma centred mesh **defineKmesh(self, NKx, Nky, Nkz)**
    - A set path **defineKmesh(self, highSymPoints)**
- Bands.x

## Figure Creation


