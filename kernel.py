import numpy as np
import auxillary as aux

from params import par
for key,val in par.items():
    exec(key + '=val')

#===============================================================================
""" M4 kernel """
#-------------------------------------------------------------------------------

def kernel_M4(s):
    w1  =   1 - (3/2)*s**2 + (3/4)**s**3
    w2  =   (1/4)*(2-s)**3
    value   =   w1*aux.heaviside(s) - w1*aux.heaviside(s,shift=1) + w2*aux.heaviside(s,shift=1) - w2*aux.heaviside(s,shift=2)
    return value/np.pi

def kernel_gradient_M4(s):
    w1  =   3*s - (9/4)*s**2
    w2  =   (3/4)*(2-s)**2
    value   =   w1*aux.heaviside(s) - w1*aux.heaviside(s,shift=1) + w2*aux.heaviside(s,shift=1) - w2*aux.heaviside(s,shift=2)
    return -value/np.pi
