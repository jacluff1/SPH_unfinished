import numpy as np
import os

params  =   {
# number of particles
'N':    1000,

# initial cloud radius
'r0':   1,

# initial rotational period
'tau0': 1,

# initla average particle speed
'vel0_av': 1e-3,

# rotational axis
'zhat': np.array([0,0,1]),

}

# def load_params():
#
#     # if not os.path.isfile("parameters.npy"): np.save("parameters.npy",params)
#     np.save("parameters.")
#     params  =   np.load("parameters.npy").item()
#     return params

# for key,val in params.items():
#     exec(key + '=val')
