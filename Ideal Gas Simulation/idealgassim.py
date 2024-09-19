import random

import numpy as np
import matplotlib.pyplot as plt
from numba import jit

class IdealGasSim():
    
    class Particle:
        ...
    
    def __init__(self, 
                 sim_size:int,
                 nparticles:int, 
                 mass:int, 
                 radius:int) -> None:
        """...

        Args:
            sim_size (int): Width and height in pixels of the simulation frame.
            nparticles (int): Starting number of particles.
            mass (int): Starting mass of particles in units of amu.
            radius (int): Starting radius of particles in pixels.
        """
        self._sim_size = sim_size
        self._nparticles = nparticles
        self._mass = mass
        self._radius = radius
        
        #TODO: Make these optional arguments.
        self._sim_depth = radius*6 # Depth in pixels of simulation.
        self._v_max = 6.8 # Maximum starting horizontal and vertical velocity of all particles in terms of pixels per timestep.
        self._vz_max = self.v_max # Maximum starting depth velocity of all particles in terms of pixels per timestep.
        self._m_per_pix = 200 # Meters per pixel.
        # Used for aggregated data collection for debugging.
        self._bin_range = 200 # Width of the velocity bins.
        self._max_range = 3000 # Max range of velocity binning.
        
        #TODO: Assign to self
        particleList = [] #the complete list of particles
        gridLength = 50 #the pixel length of each grid spacing
        gridHeight = 50 #the pixel height of each grid spacing
        length = int(simSize / gridLength + 1) #number of horizontal grid spacings
        height = int(simSize / gridHeight + 1) #number of vertical grid spacings
        grid = [[[] for _ in range(height)] for _ in range(length)] #the grid partition of particles
        
    def add_particles(n:int) -> None:
        ...
        
    def remove_particles(n:int) -> None:
        ...
        
    def time_step() -> None:
        ...
        
    def check_wall_collision(p:Particle) -> None:
        ...
        
    def wall_collision(p:Particle) -> None:
        ...
        
    def check_particle_collision(p:Particle) -> None:
        ...
        
    def particle_collision(p1:Particle, p2:Particle) -> None:
        ...
        
    def thermometer() -> float:
        ...