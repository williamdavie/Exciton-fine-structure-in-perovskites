
# Tuning the Exciton Fine Structure in 2D Perovskites

Documentation for the code supporting the completion of my Masters thesis.  

Department of Materials Science and Metallurgy,  
Cavendish Laboratory Department of Physics,  
University of Cambridge

### Contents

1. [Structure setup and relaxation](#structure-setup-and-relaxation)

## Structure setup and relaxation

### build_unit_cell.py

Sets up a unit cell with general formula Cs₂BX₄ according to distortion angles β and δ:

<div align="center">
  <img src="https://github.com/williamdavie/Exciton-fine-structure-in-perovskites/blob/main/figures/cell_marked_angles.png" width="33%">
</div>

The result satisfies the following constraints:

- Angles β and δ are fixed to a pre-set value.
- The angle between all X atoms surrounding the same B atom is 90°:
<div align="center">
  <img src="https://github.com/williamdavie/Exciton-fine-structure-in-perovskites/blob/main/figures/90degree_marked.png" width="20%">
</div>
- The bond length B–X is constant throughout the cell:

```python
BondL = np.sqrt( ( (1/4 * B_BLength**2)*(1 + np.tan(self.beta)**2) ) ) / 2
```

Once the cell is constructed it may be viewed via an output CIF file {**writeCIF(self)**} and relaxed via QE.

Relaxation details specific to our work:

- Convergence threshold of 1e-4
- 

