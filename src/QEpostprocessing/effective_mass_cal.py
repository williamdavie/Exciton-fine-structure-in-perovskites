'''

Given a band for a 2D semi-conductor, the effective mass tensor is calculated using five-point finite difference around the gamma point.

'''
import numpy as np
import math 
from read_QE_output import fetchBandData

# Example data 

filename = 'example_bandstruc.dat'

kpoints, banddata = fetchBandData(filename)
examplek = np.load('kpoints_5x5x1_gamma.npy')
exampleE = np.random.rand(len(examplek))


def effectiveMassTensor2D(kpoints: np.ndarray, Evals: np.ndarray):
    '''
    Computes the 2x2 effective mass tensor given a gamma centred 2D k point mesh

    Input: two 2D arrays containing the kpoints and corresponding energy vals
    '''

    assert kpoints.shape[0] >= 5 and kpoints.shape[1] >= 5, "5 point FD, shape must be >=5"
    
    # check gamma centred
    i_Gamma = math.floor(kpoints.shape[0]/2)
    j_Gamma = math.floor(kpoints.shape[1]/2)
    
    assert np.array(kpoints[i_Gamma,j_Gamma]).all() == np.array([0,0,0]).all(), "Must be gamma centred"

    
    
    return None

def orderMesh(kpoints: np.ndarray, Evals: np.ndarray):
    '''
    Input: a mesh of k points in shape Nx3 and corrosponding E vals in shape Nx1
    Returns: the same mesh in shape (Nx)x(Ny)x(Nz)x3 with E vals in shape (Nx)x(Ny)x(Nz)x1 
             ready for use of finite difference methods.
    '''
    
    Xvals = np.unique(kpoints[:,0])
    Yvals = np.unique(kpoints[:,1])
    Zvals = np.unique(kpoints[:,2]) # = 1 for all our cases (2D semiconductor)
    
    Nx = len(Xvals)
    Ny = len(Yvals)
    Nz = len(Zvals)
    
    kpoints3D = np.zeros((Nx,Ny,Nz,3))
    Evals3D = np.zeros((Nx,Ny,Nz))
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
            
                newKpoint = np.array([np.sort(Xvals)[i],np.sort(Yvals)[j],np.sort(Zvals)[k]])
            
                kpoints3D[i,j,k] = newKpoint 

                oldIndex = np.where(np.all(kpoints == newKpoint, axis=1))[0]
                Evals3D[i,j,k] += Evals[oldIndex][0]
    
    # reformat for 2D effective mass calculation
    
    if Nz == 1:
        return kpoints3D.reshape(Nx,Ny,3),Evals3D.reshape(Nx,Ny)
    else:
        return kpoints3D,Evals3D
        
        
kvals,Evals = orderMesh(examplek,exampleE)

effectiveMassTensor2D(kvals,Evals)


    
    