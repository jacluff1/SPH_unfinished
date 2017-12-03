import djak.phys.SPH.kernel as kernel
import djak.math as dm
import numpy as np

#===============================================================================
""" parameters """
#-------------------------------------------------------------------------------

c0      =   1
# rho0    =   1e-14
rho0    =   10000
alpha   =   1
beta    =   1
gamma   =   1
G       =   1e-4

#===============================================================================
""" sound speed in gas """
#-------------------------------------------------------------------------------

def c_gas(j,rhoA,PA):
    """ compute isothermal sound speed in particle (and approximates adiabadic)

    args
    ----
    j:      index of particle j
    rhoA:   array of particle densities
    PA:     array of particle pressures
    """
    return np.sqrt( PA[j] / rhoA[j] )

#===============================================================================
""" density """
#-------------------------------------------------------------------------------

def rho(j,rA,mA,hA,W=kernel.W_M4):
    """ computes density of particle i

    args
    ----
    j:  particle index
    rA: array of particle positions
    mA: array of particle masses
    hA: array of smoothing lengths
    W:  ** smoothing function - default = M4
    """
    assert rA.shape[0] == mA.shape[0] == hA.shape[0], "arrays are not matched"
    N   =   len(mA)

    tot =   0
    for i in range(N):
        if i != j:

            r_ij1   =   np.linalg.norm( rA[j,:] - rA[i,:] )
            m_i     =   mA[i]
            h_ij    =   .5 * (hA[i] + hA[j])

            tot     +=  m_i * h_ij**(-3) * W(r_ij1,h_ij)

    return tot

#===============================================================================
""" presure """
#-------------------------------------------------------------------------------

def P_bol(j,rhoA):
    """ compute barotropic pressure from density

    args
    ----
    j:      index of particle
    rhoA:   array of densities
    """

    rho_j   =   rhoA[j]
    return rho_j * c0**2 * np.sqrt( 1 + (rho_j/rho0)**(4/3) )

#===============================================================================
""" acelleration terms """
#-------------------------------------------------------------------------------

def acc_fluid(j,rA,mA,rhoA,PA,hA,dW=kernel.dW_M4):
    """ compute acceleration from fluid dynamics

    args
    ----
    j:      particle index
    rA:     array of particle positions
    mA:     array of particle masses
    rhoA:   array of particle densities
    PA:     array of particle pressures
    hA:     array of particle smoothing lengths
    dW:     ** smoothing function - default = kernel.dW_M4
    """
    assert rA.shape[0] == mA.shape[0] == rhoA.shape[0] == PA.shape[0] == hA.shape[0], "arrays are mismatched"
    N       =   len(mA)
    rho_j   =   rhoA[j]
    P_j     =   PA[j]

    tot     =   0
    for i in range(N):
        if i != j:

            r_ij    =   rA[j,:] - rA[i,:]
            r_ij1   =   np.linalg.norm(r_ij)
            m_i     =   mA[i]
            rho_i   =   rhoA[i]
            P_i     =   PA[i]
            h_ij    =   0.5 * (hA[i] + hA[j])

            # if np.array_equal(r_ij,np.zeros(3)): print(i,j)
            # assert r_ij1 != 0, "r_ij1"
            # assert h_ij != 0, "h_ij"
            # assert rho_i != 0, "rho_i"
            # assert rho_j != 0, "rho_j"

            tot     +=  m_i * (r_ij/r_ij1) * h_ij**(-4) * ( (P_i/rho_i**2) + (P_j/rho_j**2) ) * dW(r_ij1,h_ij)

    return - tot

def acc_visc(j,rA,vA,mA,rhoA,PA,hA,dW=kernel.dW_M4):
    """ compute acceleration from artificial viscosity

    args
    ----
    j:  particle index
    rA: array of particle positions
    vA: array of particle velocities
    hA: array of smoothing lengths
    dW: ** smoothing function - default = kernel.dW_M4
    """
    assert rA.shape[0] == vA.shape[0] == mA.shape[0] == rhoA.shape[0] == hA.shape[0], "arrays are not matched"
    N       =   len(mA)
    c_j     =   c_gas(j,rhoA,PA)

    tot     =   0
    for i in range(N):
        if i != j:

            r_ij    =   rA[j,:] - rA[i,:]
            r_ij1   =   np.linalg.norm(r_ij)
            v_ij    =   vA[j,:] - vA[i,:]
            m_i     =   mA[i]
            c_i     =   c_gas(i,rhoA,PA)
            c_ij    =   0.5 * (c_i + c_j)
            h_ij    =   0.5 * (hA[i] + hA[j])
            rho_ij  =   0.5 * (rhoA[i] + rhoA[j])

            c       =   np.dot(v_ij,r_ij)
            mu_ij   =   ( c * h_ij ) / ( r_ij1**2 + 0.01*h_ij**2 )

            a       =   ( -alpha * mu_ij * c_ij + beta * mu_ij**2 ) / rho_ij
            b       =   0
            Pi_ij   =   a*dm.heavi(-c) + b*dm.heavi(c)

            # if Pi_ij == 0:
            #     print("i,j:",i,j)
            #     print("c:",c)
            #     print("c_ij",c_ij)
            #     print("")
            #     assert Pi_ij != 0

            tot     +=  m_i * h_ij**(-4) * Pi_ij * dW(r_ij1,h_ij) * (r_ij/r_ij1)

    return - tot

def acc_grav(j,rA,mA,hA,W=kernel.W_M4star):
    """ compute acceleration on particle j from gravity

    args
    ----
    j:  index of particle j
    rA: array of particle positions
    mA: array of particle masses
    hA: array of particle smoothing lengths
    W:  ** smoothing function - default = kernel.W_M4star
    """
    assert rA.shape[0] == mA.shape[0] == hA.shape[0], "arrays are mismatched"
    N       =   len(mA)

    tot     =   0
    for i in range(N):
        if i != j:

            r_ij    =   rA[j,:] - rA[i,:]
            r_ij1   =   np.linalg.norm(r_ij)
            m_i     =   mA[i]
            h_ij    =   0.5 * (hA[i] + hA[j])

            tot     +=  m_i * (r_ij/r_ij**3) * W(r_ij1,h_ij)

    return - tot * G

def acc_total(j,rA,vA,mA,rhoA,PA,hA,Wstar=kernel.W_M4star,dW=kernel.dW_M4):

    # a_fluid =   acc_fluid(j,rA,mA,rhoA,PA,hA,dW=dW)
    # a_visc  =   acc_visc(j,rA,vA,mA,rhoA,PA,hA,dW=dW)
    a_grav  =   acc_grav(j,rA,mA,hA,W=Wstar)

    # return a_fluid + a_visc + a_grav
    # return a_grav
