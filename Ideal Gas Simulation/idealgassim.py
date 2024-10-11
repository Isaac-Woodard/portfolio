import random
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
    
@dataclass
class Particle:
    P: np.ndarray = np.zeros(3)
    V: np.ndarray = np.zeros(3)
    m: float
    r: float
    canvasID = None #TODO: Still useful after refactor?
    tagged = None #TODO: Still useful after refactor?
        
    @property
    def speed(self):
        return np.sum(np.power(self.V, 2))**0.5
    
    @property
    def KE(self):
        s = self.speed
        return 0.5*self.m*s**2
    
    @property
    def momentum(self):
        s = self.speed
        return self.m*s


class IdealGasSim():
    def __init__(self, 
                 sim_size:int,
                 nparticles:int, 
                 mass:int, 
                 radius:int,
                 external_temp:float,
                 vmax:float=6.8,
                 m_per_pix:float=200,
                 vbin_width:float=200,
                 vbin_range:float=3000,
                 grid_size:int=50) -> None:
        """Creates a 3D ideal gas in a box. The gas is modeled as a group of discrete particles obeying classical mechanics.

        Args:
            sim_size (int): Width and height in pixels of the simulation frame.
            nparticles (int): Starting number of particles.
            mass (int): Starting mass of particles in units of amu.
            radius (int): Starting radius of particles in pixels.
            external_temp (float): Temperature of environment in Kelvin.
            vmax (float): Max starting velocity of particles in pixels per timestep.
            m_per_pix (float): Meters per pixel.
            vbin_width (float): Width of velocity bins used for data aggregation.
            vbin_range (float): Range for velocity binning.
            grid_size (int): Pixel length of particle list partitions.
        """
        self._sim_size = sim_size
        self._nparticles = nparticles
        self._mass = mass
        self._radius = radius
        self._external_temp = external_temp
        
        self._sim_depth = radius*6 # Depth in pixels of simulation.
        self._vmax = vmax
        self._vz_max = vmax # Maximum starting depth velocity.
        self._m_per_pix = m_per_pix
        self._vbin_width = vbin_width
        self._vbin_range = vbin_range
        
        self.particles = [] #TODO: Add some starting particles?
        self._grid_size = grid_size
        n = int(sim_size / grid_size + 1)
        self.partition = [[[] for _ in range(n)] for _ in range(n)]
        
    def add_particles(self, n:int) -> None:
        """Adds particles to the simulation.

        Args:
            n (int): Number of particles to add.
        """
        for _ in range(n):
            dx = random.uniform(-self._vmax, self._vmax)
            dy = random.uniform(-self._vmax, self._vmax)
            dz = random.uniform(-self._vmax, self._vmax)
            overlap = True
            while overlap: # Check for overlap with existing particles.
                x = random.randint(1+p.r, self._sim_size-p.r)
                y = random.randint(1+p.r, self._sim_size-p.r)
                z = random.randint(1+p.r, self._sim_depth-p.r)
                for p in self.particles:
                    if x >= p.p[0]-2*p.r and x <= p.p[0]+2*p.r and y >= p.p[1]-2*p.r and y <= p.p[1]+2*p.r and z >= p.p[2]-2*p.r and z <= p.p[2]+2*p.r:
                        overlap = True
                        break
                else:
                    overlap = False 
                    
            p = Particle(P=np.ndarray((x, y, z)), 
                         V=np.ndarray((dx, dy, dz)),
                         m=self._mass, r=self._radius)    
            nx = int(p.P[0] / self._grid_size)
            ny = int(p.P[1] / self._grid_size)
            self.partition[nx][ny].append(p)
            self.particles.append(p)
        
    def remove_particles(self, n:int) -> None:
        """Removes particles from the simulation.

        Args:
            n (int): Number of particles to remove.
        """
        for _ in range(n):
            self.particles.pop()
        
    #TODO: Refactor
    def time_step(self, n:int=1) -> None:
        """Advances the simulation n time-steps.
        #TODO Describe the units for the time-steps.

        Args:
            n (int, optional): Number of time steps. Defaults to 1.
        """
        if unpaused is True:
            #clears the grid
            i = int(self._sim_size / self._grid_size + 1)
            self.partition = [[[] for _ in range(n)] for _ in range(n)]
            #update each particle position and the grid
            dt = simSpeed.get() / 50 #parameter to adjust simulation speed
            for p in particleList:
                p.p[0] = p.p[0] + p.v[0] * dt
                p.p[1] = p.p[1] + p.v[1] * dt
                p.p[2] = p.p[2] + p.v[2] * dt #3D Edit Marker
                gridx = int(p.p[0] / gridLength)
                gridy = int(p.p[1] / gridHeight)
                grid[gridx][gridy].append(p)
            #check for collisions with wall and particles and update trajectories if necessary
            for p in particleList:
                if checkWallCollision(p):
                    wallCollision(p)
                else:
                    p2 = checkParticleCollision(p) #store the colliding particle
                    if p2 != p:
                        particleCollision(p, p2)
        
        #TODO: This should probably be its own function.
        #code to automate bin counts
        global binTimes
        binNum = 0
        for i in range(0, maxRange, binRange):
            count = 0
            #get high and low speed values from text fields
            low = i
            high = i+binRange
            #step through particles and find speeds
            for p in particleList:
                speed = mPerPixel * (p.v[0]**2 + p.v[1]**2 +p.v[2]**2)**0.5 #normalizing factor: 200 meters per pixel
                if speed >= low and speed <= high:
                    count += 1
            binCounts[binNum] = (binCounts[binNum]*binTimes+count)/(binTimes+1)
            binNum += 1
        binTimes += 1
        
    #TODO: Return which wall the collision is with. Use info for wall_collision() call.
    def check_wall_collision(self, p:Particle) -> bool:
        """Checks the particle for wall collisions.

        Args:
            p (Particle): The particle to check.

        Returns:
            Bool: Whether there is a wall collision.
        """
        if p.P[0]-p.r < 1 or p.P[0]+p.r > self._sim_size:
            return True
        if p.P[1]-p.r < 1 or p.P[1]+p.r > self._sim_size:
            return True
        if p.P[2]-p.r < 1 or p.P[2]+p.r > self._sim_depth:
            return True
        else:
            return False
        
    def wall_collision(self, p:Particle) -> None:
        """Updates particle velocity in case of a wall collision.

        Args:
            p (Particle): The particle to update.
        """
        # Bounce the particle.
        if p.P[0]-p.r < 1 or p.P[0]+p.r > self._sim_size: 
            p.V[0] = -p.V[0]
            # Make sure the particle isn't pushed behind the wall.
            if p.P[0]-p.r < 1:
                p.P[0] = p.r + 1
            if p.P[0]+p.r > self._sim_size:
                p.P[0] = self._sim_size - 2 - p.r
        if p.P[1]-p.r < 1 or p.P[1]+p.r > self._sim_size:
            p.V[1] = -p.V[1]
            if p.P[1]-p.r < 1:
                p.P[1] = p.r + 1
            if p.P[1]+p.r > self._sim_size:
                p.P[1] = self._sim_size - 2 - p.r
        if p.P[2]-p.r < 1 or p.P[2]+p.r > self._sim_depth:
            p.V[2] = -p.V[2]
            if p.P[2]-p.r < 1:
                p.P[2] = p.r + 1
            if p.P[2]+p.r > self._sim_depth:
                p.P[2] = self._sim_depth - 2 - p.r
            
        # Simulate heating or cooling the particle to meet target temperature.
        k = 1000
        if temperature != self._external_temp: #NOTE: Might be more efficient to skip the check.
            p.V = p.V + (self._external_temp - temperature) / k * p.V
        
    def check_particle_collision(self, p:Particle) -> Particle:
        """Checks for particle collisions. Search is confined to neighboring partitions.
        Does not check for collisions with multiple particles.

        Args:
            p (Particle): Particle to check.

        Returns:
            Particle: The colliding particle if found, otherwise None.
        """
        gridx = int(p.P[0] / self._grid_size)
        gridy = int(p.P[1] / self._grid_size)
        # Step through adjacent grid spaces as well.
        for i in (-1,0,1):
            for j in (-1,0,1):
                if gridx+i < 0 or gridy+j < 0 or gridx+i >= self._sim_size or gridy+j >= self._sim_size:
                    continue
                for p2 in grid[gridx+i][gridy+j]:
                    if p is p2:
                        continue
                    dmin = p.r + p2.r
                    D = p.P - p2.P
                    if (abs(D[0]) > dmin or abs(D[1]) > dmin or abs(D[2]) > dmin):
                        continue
                    d = (D[0]**2 + D[1]**2 + D[2]**2)**0.5
                    if (d > dmin):
                        continue
                    return p2
        return None
        
    def particle_collision(self, p1:Particle, p2:Particle) -> None:
        """Update velocities of colliding particles p1 and p2.

        Args:
            p1 (Particle): The first particle.
            p2 (Particle): The second particle.
        """
        # Algorithm modified from Gas.java (C) 2001 by Paul Falstad, www.falstad.com
        # Go to http://exploratoria.github.io/exhibits/mechanics/elastic-collisions-in-3d/ for more details on the theory
     
        # Calculate time since collision with quadratic formula: 
        # [(x1-x2)+t(dx1-dx2)]^2 + [(y1-y2)+...]^2 + [(z1-z2)+...]^2 = dmin^2
        Pdif = p1.P - p2.P
        Vdif = p1.V - p2.V
        dmin = p1.r + p2.r
        a = np.sum(np.power(Vdif, 2))
        b = 2 * np.sum(np.multiply(Pdif, Vdif))
        c = np.sum(np.power(Pdif, 2)) - dmin**2
        t = (-b - (b**2 - 4*a*c)**0.5) / (2*a) # Roots
        t2 = (-b + (b**2 - 4*a*c)**0.5) / (2*a)
        if (abs(t) > abs(t2)):
            t = t2
        
        # Return p1 to place and time of collision
        p1.P += t * p1.V

        # Find the unit vector components between the two particles' centers to treat the collision as one-dimensional.
        V = p1.P - p2.P
        norm = np.sqrt(np.sum(np.power(V, 2)))
        Vnorm = V/norm

        # Treat the particles as one particle and calculate the center of mass velocity.
        totmass = p1.m + p2.m
        COMV = (p1.m * p1.V + p2.m * p2.V) / totmass

        # Find momentum transferred to p1.
        pn = np.multiply(Vnorm, p1.V-COMV)
        P = 2 * Vnorm * pn

        # Transfer momentum to p1.
        p1.V -= P

        # Remove momentum from p2.
        massratio = p1.m/p2.m
        p2.V += P * massratio
        
    def thermometer(self) -> float:
        """Calculate the temperature of the gas based on the Maxwell Boltzmann distribution.

        Returns:
            float: The temperature in Kelvin.
        """
        KEtotal = 0
        for p in self.particles:
            speed = self._m_per_pix * p.speed # 200 meters per pixel
            KEtotal += 0.5 * p.m * 1.66e-27 * speed**2 # proton mass 1.67e-27 kg
        KEmean = KEtotal / len(self.particles)
        temp = KEmean / 1.38e-23 / 1.5
        return temp
    
    #TODO: Plot current sim state.
    def plot(self, figure:plt.Figure=None) -> None:
        ...