import numpy as np

from params import par
for key,val in par.items():
    exec(key + '=val')

def soundspeed_isothermal(T):
    return np.sqrt( k * T / mu )

def pressure(soundspeed,density):
    return density * soundspeed**2 * np.sqrt(1 + (density/rho0)**(4/3) )
