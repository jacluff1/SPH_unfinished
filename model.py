import numpy as np
import initialConditions as ic
import SPH as sph
import physicalEquations as pe
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as ani

from params import par
for key,val in par.items():
    exec(key + '=val')

# class particle:
#     def __init__(self,position,velocity,mass,smoothingLength):
#         self.r  =   position
#         self.v  =   velocity
#         self.m  =   mass
#         self.h  =   smoothingLength

def run_particle_motion():
    rmax    =   r0*1.2

    # arrays that aren't updated
    Mass            =   np.ones(N) * m

    # position and velocity arrays
    Position        =   ic.placement_uniform_sphere()
    Velocity        =   ic.velocity_zero()
    # Velocity        =   ic.velocity_random_motion() + ic.velocity_rotation_uniform(Position)
    # Velocity        =   ic.velocity_rotation_uniform(Position)

    # use for now, but requires further development
    dt              =   .001

    FFMpegWriter    =   ani.writers['ffmpeg']
    metadata        =   dict(title='SPH Star Formation', artist='Matplotlib')
    writer          =   FFMpegWriter(fps=15, metadata=metadata)

    fig =   plt.figure(figsize=[15,15])
    ax  =   fig.gca(projection='3d')
    ax.set_aspect(1)
    ax.plot(Position[:,0],Position[:,1],Position[:,2], 'go')

    with writer.saving(fig,"SPH.mp4", 100):
        for i in range(100):

            ax.clear()

            # use RK2
            Acceleration    =   sph.acceleration_total(Position,Velocity,Mass)
            Position_half   =   Position + Velocity*(dt/2)
            Velocity_half   =   Velocity + Acceleration*(dt/2)

            Acceleration_half   =   sph.acceleration_total(Position_half,Velocity_half,Mass)
            Position            +=  Velocity_half*dt
            Velocity            +=  Acceleration_half*dt

            ax.plot(Position[:,0],Position[:,1],Position[:,2], 'go')
            ax.set_aspect(1)
            ax.set_xlim([-rmax,rmax])
            ax.set_ylim([-rmax,rmax])
            ax.set_zlim([-rmax,rmax])
            writer.grab_frame()

    plt.close()
