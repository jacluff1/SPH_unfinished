import djak.math as dm
import numpy as np

heavi   =   dm.heavi

#===============================================================================
""" M4 kernel """
#-------------------------------------------------------------------------------

def W_M4(r_mag,h):
    """ M4 smoothing kernel

    args
    ----
    r_mag:  magnitude of position vector
    h:      smoothing length
    """

    s   =   r_mag/h

    a   =   1 - (3/2)*(s**2) + (3/4)*(s**3)
    b   =   (1/4)*(2 - s)**3
    c   =   0

    return (1/np.pi) * ( a*heavi(s)*heavi(1-s) + b*heavi(s-1)*heavi(2-s) + c*heavi(s-2) )

def dW_M4(r_mag,h):
    """ derivative of M4 smoothing kernel

    args
    ----
    r_mag:  magnitude of position vector
    h:      smoothing length
    """

    s   =   r_mag/h

    a   =   3*s - (9/4)*s**2
    b   =   (3/4)*(2-s)**2
    c   =   0

    return -(1/np.pi) * ( a*heavi(s)*heavi(1-s) + b*heavi(s-1)*heavi(2-s) + c*heavi(s-2) )

def W_M4star(r_mag,h):
    """ M4 smoothing kernel for Tree Code Gravity

    args
    ----
    r_mag:  magnitude of position vector
    h:      smoothing length
    """

    s   =   r_mag/h

    a   =   40*s**3 - 36*s**5 + 15*s**6
    b   =   80*s**3 - 90*s**4 + 36*s**5 - 5*s**6 - 2
    c   =   30

    return (1/30) * ( a*heavi(s)*heavi(1-s) + b*heavi(s-1)*heavi(2-s) + c*heavi(s-2) )
