'''

Given a band for a 2D semi-conductor, the effective mass tensor is calculated using five-point finite difference around the gamma point.

'''
import numpy as np
import math 
from read_QE_output import *

# Example data 


def effectiveMassTensor2D(kpoints: np.ndarray, Evals: np.ndarray, celldims: np.ndarray):
    '''
    Computes the 2x2 effective mass tensor given a gamma centred 2D k point mesh

    Input: two 2D arrays containing the kpoints and corresponding energy vals
    
            + the celldims to correctly scale dk
    
    '''

    assert kpoints.shape[0] >= 5 and kpoints.shape[1] >= 5, "5 point FD, shape must be >=5"
    
    # check gamma centred
    i_Gamma = math.floor(kpoints.shape[0]/2)
    j_Gamma = math.floor(kpoints.shape[1]/2)
    
    assert np.array(kpoints[i_Gamma,j_Gamma]).all() == np.array([0,0,0]).all(), "Must be gamma centred"
    
    # use formula for 5 point difference
    
    # assumes evenly spaced
    dkx = round((kpoints[1,0] - kpoints[0,0])[0],10)*(2*np.pi/celldims[0])
    dky = round((kpoints[0,1] - kpoints[0,0])[1],10)*(2*np.pi/celldims[1])
    
    
    #------------- d2/dx2 -------------
    Exvals = Evals[i_Gamma,j_Gamma-2:j_Gamma+3].reshape(5)
    
    d2dx2 = fivePoint2ndDeriv(Exvals,dkx)

    #------------- d2/dy2 -------------
    Eyvals = Evals[i_Gamma-2:i_Gamma+3,j_Gamma].reshape(5)
    
    d2dy2 = fivePoint2ndDeriv(Eyvals,dky)
    
    #-------- d2/dxdy = d2/dydx -------
    Exymesh = Evals[i_Gamma-2:i_Gamma+3,j_Gamma-2:j_Gamma+3].reshape(5,5)
    
    d2dxdy = d2dydx = MixedDeriv(Exymesh,dkx,dky) 
    
    M = np.zeros((2,2))
    
    M[0,0] = d2dx2
    M[1,1] = d2dy2
    M[0,1] = d2dxdy
    M[1,0] = d2dydx
    
    M = 1/(6.582e-16)**2*M*((1e-10)**2)*(1/(1.602e-19))
    
    return M

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
    
def fivePoint2ndDeriv(points: np.ndarray, stepsize: float) -> float:
    '''
    Uses a central difference to calculate the second derivative at the centre of 5 points
    '''
    assert points.shape[0]==5, 'must input 1D array of 5 points'
    F = np.copy(points) # formatting
    result = 1/(12 * stepsize**2) * ( - (F[0] + F[4]) 
                                     + 16 * ( F[1] + F[3]) 
                                     - 30 * F[2])
    
    return result


def MixedDeriv(points: np.ndarray, stepsize1: float, stepsize2: float) -> float:
    '''
    Uses a central difference to calculate the mixed derivative
    '''
    assert points.shape[0:2] == (5,5), 'must input a 5x5 array'
    F = np.copy(points) 
    result = 1/(4*stepsize1*stepsize2) * (F[3,3] - F[3,1] - F[1,3] + F[1,1])
    return result
    