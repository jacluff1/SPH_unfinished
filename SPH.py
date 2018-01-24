import numpy as np
import physicalEquations as pe
import auxillary as aux

import kernel as kernel
W   =   kernel.kernel_M4
Wg  =   kernel.kernel_gradient_M4

from params import par
for key,val in par.items():
    exec(key + '=val')

#===============================================================================
""" density """
#===============================================================================

def density_j(j,R,M,H):
    total   =   0
    for i in range(N):
        r_ij    =   R[j,:] - R[i,:]
        h_ij    =   .5 * (H[j] + H[i])
        s_ij    =   np.linalg.norm(r_ij)/h_ij
        total   +=  M[i] * h_ij**(-3) * W(s_ij)
    return total

def density_j_dot(j,R,V,M,H,P):
    total   =   0
    for i in range(N):

        r_ij    =   R[j,:] - R[i,:]
        v_ij    =   V[j,:] - V[i,:]
        h_ij    =   .5 * (H[j] + H[i])
        s_ij    =   np.linalg.norm(r_ij)/h_ij

        vec1    =   M[i] * h_ij**(-4) * v_ij
        vec2    =   Wg(s_ij) * aux.direction(r_ij)

        total   +=  np.dot( vec1 , vec2 )

#===============================================================================
""" artifitial viscosity """
#===============================================================================

def acceleration_visc_j(j,R,V,M,H,Rho,C):
    total   =   0
    for i in range(N):
        if i != j:
            r_ij    =   R[j,:] - R[i,:]
            v_ij    =   V[j,:] - V[i,:]
            h_ij    =   .5 * (H[j] + H[i])
            s_ij    =   np.linalg.norm(r_ij)/h_ij
            c_ij    =   .5 * (C[j] + C[i])
            mu_ij   =   (np.dot(v_ij,r_ij)*h_ij) / (np.linalg.norm(r_ij)**2 + 0.01*h_ij**2)
            rho_ij  =   .5 * (Rho[j] + Rho[i])

            if v_ij.dot(r_ij) < 0:
                Pi_ij   =   (-alpha*mu_ij*c_ij + beta*mu_ij**2) / rho_ij
            else:
                Pi_ij   =   0

            total   +=  M[i] * h_ij**(-4) * Pi_ij * Wg(s_ij) * aux.direction(r_ij)
    return -total

#===============================================================================
""" temporary gravity, will replace with "tree-code" system in tree.py """
#===============================================================================

def acceleration_grav_j(j,R,M):
    total   =   0
    for i in range(N):
        if i != j:
            r_ij    =   R[j,:] - R[i,:]
            total   +=  M[i]/np.linalg.norm(r_ij)**2 * aux.direction(r_ij)
    return -G*total

#===============================================================================
""" particle acceleration """
#===============================================================================

def velocity_j_dot(j,R,M,H,Rho,P,Agrav,Avisc):
    total   =   0
    for i in range(N):
        if i != j:
            r_ij    =   R[j,:] - R[i,:]
            h_ij    =   .5 * (H[j] + H[i])
            s_ij    =   np.linalg.norm(r_ij)/h_ij
            total   +=  M[i] * h_ij**(-4) * ( (P[i]/Rho[i]**2) + (P[j]/Rho[j]**2) ) * Wg(s_ij) * aux.direction(r_ij)
    return -total + Agrav[j] + Avisc[j]

def acceleration_total(Position,Velocity,Mass):

    # not implemented, use for now
    # Temperature         =   np.ones(N) * 1
    SoundSpeed          =   np.ones(N) * 1e-1
    SmoothingLength     =   np.ones(N) * r0

    # SmoothingLength - not implemented
    Density             =   np.array([ density_j(j,Position,Mass,SmoothingLength) for j in range(N) ])
    # Temperature - not implemented
    # SoundSpeed - not implemented
    Pressure            =   pe.pressure(SoundSpeed,Density)
    Acceleration_visc   =   np.array([ acceleration_visc_j(j,Position,Velocity,Mass,SmoothingLength,Density,SoundSpeed) for j in range(N) ])
    Acceleration_grav   =   np.array([ acceleration_grav_j(j,Position,Mass) for j in range(N) ])

    return np.array([ velocity_j_dot(j,Position,Mass,SmoothingLength,Density,Pressure,Acceleration_grav,Acceleration_visc) for j in range(N) ])
