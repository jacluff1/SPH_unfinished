import numpy as np

from params import params
for key,val in params.items():
    exec(key + '=val')

def placement_uniform_sphere():
    """
    make a random distribution of 'N' particles
    inside a sphere.

    returns
    -------
    2D position array

    thanks
    ------
    https://stackoverflow.com/a/5408843
    """

    # set up 1D position arrays in SPC
    U           =   np.random.uniform(0,1,N)
    COSTHETA    =   np.random.uniform(-1,1,N)

    R       =   r0 * U**(1/3)
    THETA   =   np.arccos(COSTHETA)
    PHI     =   np.random.uniform(0,2*np.pi,N)

    # set up 2D position array
    X       =   np.zeros( (N,3) )
    # convert SPC to CC
    X[:,0]  =   R * np.sin(THETA) * np.cos(PHI)
    X[:,1]  =   R * np.sin(THETA) * np.sin(PHI)
    X[:,2]  =   R * np.cos(THETA)

    return X

def velocity_random_motion():

    V1d =   vel0_av * np.random.randn(N)

    return V1d
