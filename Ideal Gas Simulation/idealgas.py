from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from numba import jit

from numpy import ndarray
    
@dataclass
class Particle:
    P: ndarray # Position in pixels
    V: ndarray # Velocity in pixels per timestep
    m: float # Mass in amu
    r: float # Radius in pixels
        
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
        """Creates a 3D ideal gas in a box. The gas is modeled as a group of 
        discrete particles obeying classical mechanics.

        Args:
            sim_size (int): Width and height in pixels of the simulation frame.
            nparticles (int): Starting number of particles.
            mass (int): Starting mass of particles in units of amu.
            radius (int): Starting radius of particles in pixels.
            external_temp (float): Temperature of environment in Kelvin.
            vmax (float): Max starting velocity of particles along each axis in pixels per timestep.
            m_per_pix (float): Meters per pixel.
            TODO: Might need to impose a dt for a timestep.
            vbin_width (float): Width of velocity bins used for data aggregation.
            vbin_range (float): Range for velocity binning.
            grid_size (int): Pixel length of particle list partitions.
        """
        self._sim_size = np.asarray((sim_size, sim_size, radius*6)) # X Y Z
        self._nparticles = nparticles
        self._mass = mass
        self._radius = radius
        self._external_temp = external_temp
        
        self._Vmax = np.asarray((vmax, vmax, vmax))
        self._m_per_pix = m_per_pix
        self._vbin_width = vbin_width
        self._vbin_range = vbin_range
        
        self.particles = []
        self._grid_size = grid_size
        n = int(sim_size / grid_size + 1)
        self.partition = [[[] for _ in range(n)] for _ in range(n)]
        self._i_step = 0
        
        if nparticles:
            self.add_particles(nparticles)
        
    def add_particles(self, n:int) -> None:
        """Adds particles to the simulation.

        Args:
            n (int): Number of particles to add.
        """
        for _ in range(n):
            V = np.random.uniform(-self._Vmax, self._Vmax)
            overlap = True
            i = 0
            while overlap and i < 100: # Check for overlap with existing particles.
                P = np.random.randint(1+self._radius, self._sim_size-self._radius)
                for p in self.particles:
                    if np.any(P >= p.P-2*p.r) and np.any(P <= p.P+2*p.r):
                        overlap = True
                        i +=1
                        break
                overlap = False
                    
            p = Particle(P=P, V=V, m=self._mass, r=self._radius)    
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
        
    def time_step(self, n:int=1) -> None:
        """Advances the simulation n time-steps.

        Args:
            n (int, optional): Number of time steps. Defaults to 1.
        """
        for _ in range(n):
            # Reset the partition.
            i = int(self._sim_size[0] / self._grid_size + 1)
            self.partition = [[[] for _ in range(i)] for _ in range(i)]
            
            # Update particle positions.
            for p in self.particles:
                p.P = p.P + p.V
                nx = int(p.P[0] / self._grid_size)
                ny = int(p.P[1] / self._grid_size)
                self.partition[nx][ny].append(p)
                
            # Check for collisions.
            for p in self.particles:
                if self.check_wall_collision(p):
                    self.wall_collision(p)
                else:
                    p2 = self.check_particle_collision(p)
                    if p2:
                        self.particle_collision(p, p2)
        self._i_step += n
        
    def speed_bin(self, min_speed:float, max_speed:float) -> int:
        """Return the number of particles within a speed range.

        Args:
            min_speed (float): The minimum speed, inclusive.
            max_speed (float): The maximum speed, exclusive.

        Returns:
            int: The number of particles.
        """
        count = 0
        for p in self.particles:
            speed = self._m_per_pix * (p.V[0]**2 + p.V[1]**2 +p.V[2]**2)**0.5
            if speed >= min_speed and speed < max_speed:
                count += 1
        return count
        
    def check_wall_collision(self, p:Particle) -> np.ndarray:
        """Checks the particle for wall collisions.

        Args:
            p (Particle): The particle to check.

        Returns:
            ndarray: Whether there is a wall collision along the X, Y, and Z axes.
        """
        return p.P-p.r < 1 | p.P+p.r > self._sim_size
        
    def wall_collision(self, p:Particle, checks:np.ndarray) -> None:
        """Updates particle velocity in case of a wall collision.
        
        NOTE: The particle postion isn't backed up the 
        way it is in a particle collision, so it will effecticely 
        pass the wall before bouncing.

        Args:
            p (Particle): The particle to update.
            checks (ndarray): An array indicating whether there was a wall collision along the X, Y, and Z axes.
        """
        # Bounce the particle.
        i = np.where(checks == 1)
        if checks[i]: 
            p.V[i] = -p.V[i]
            # Make sure the particle isn't pushed behind the wall.
            if p.P[i]-p.r < 1:
                p.P[i] = p.r + 1
            if p.P[i]+p.r > self._sim_size[i]:
                p.P[i] = self._sim_size[i] - 2 - p.r
            
        # Simulate heating or cooling the particle to meet target temperature.
        k = 1000 # Arbitrary constant.
        p.V = p.V + (self._external_temp - self.thermometer() ) * p.V / k 
        
    def check_particle_collision(self, p:Particle) -> Particle:
        """Checks for particle collisions. Search is confined to neighboring partitions.
        Does not check for collisions with multiple particles.

        Args:
            p (Particle): Particle to check.

        Returns:
            Particle: The colliding particle if found, otherwise None.
        """
        nx = int(p.P[0] / self._grid_size)
        ny = int(p.P[1] / self._grid_size)
        # Step through adjacent grid spaces as well.
        for i in (-1,0,1):
            for j in (-1,0,1):
                if nx+i < 0 or ny+j < 0 or nx+i >= self._sim_size[0] or ny+j >= self._sim_size[0]:
                    continue
                for p2 in self.partition[nx+i][ny+j]:
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
        #TODO: Describe algorithm in a separate markdown file.
        # Calculate time since collision with quadratic formula: 
        # [(x1-x2)+t(dx1-dx2)]^2 + [(y1-y2)+t(dy1-dy2)]^2 + [(z1-z2)+t(dz1-dz2)]^2 = d_min^2
        # Ref: Algorithm modified from Gas.java (C) 2001 by Paul Falstad, www.falstad.com
        P_dif = p1.P - p2.P
        V_dif = p1.V - p2.V
        d_min = p1.r + p2.r
        a = np.sum(np.power(V_dif, 2))
        b = 2 * np.sum(np.multiply(P_dif, V_dif))
        c = np.sum(np.power(P_dif, 2)) - d_min**2
        t = (-b - (b**2 - 4*a*c)**0.5) / (2*a)
        t2 = (-b + (b**2 - 4*a*c)**0.5) / (2*a)
        if (abs(t) > abs(t2)):
            t = t2
        
        # Return p1 to place and time of collision
        p1.P += t * p1.V

        # Given momentum and kinetic energy are conserved, solve for final velocities.
        # Ref: https://en.wikipedia.org/wiki/Elastic_collision
        p1.V = p1.V*((p1.m - p2.m)/(p1.m + p2.m)) + p2.V*(2*p2.m/(p1.m + p2.m))
        p2.V = p1.V*(2*p1.m/(p1.m + p2.m)) + p2.V*((p2.m - p1.m)/(p1.m + p2.m))
        
        # # Find the unit vector components between the two particles' centers to treat the collision as one-dimensional.
        # V = p1.P - p2.P
        # norm = np.sqrt(np.sum(np.power(V, 2)))
        # V_norm = V/norm

        # # Treat the particles as one particle and calculate the center of mass velocity.
        # V_COM = (p1.m * p1.V + p2.m * p2.V) / (p1.m + p2.m)

        # # Find momentum transferred to p1.
        # pn = np.multiply(V_norm, p1.V-V_COM)
        # V_P = 2 * V_norm * pn

        # # Transfer momentum.
        # p1.V -= V_P * p2.m/p1.m
        # p2.V += V_P * p1.m/p2.m
        
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
    
    def plot(self, figure:plt.Figure=None) -> None:
        """Plot the current positions of the particles.
        
        #NOTE: Could easily add a 3d option with matplotlib's 3d scatter plot.

        Args:
            figure (plt.Figure, optional): Figure to plot on. Will create a new figure if None.
        """
        show = False
        if figure is None:
            figure = plt.figure()
            show = True
        
        ax = figure.add_subplot()
        
        for p in self.particles:
            radius = p.r * p.P[2]/self._sim_size[2]
            circle = plt.Circle((p.P[0], p.P[1]), radius, color="b")
            ax.add_artist(circle)
            
        ax.set(title=f"Ideal Gas Simulation: Timestep {self._i_step}", 
               xlabel="X Pos", 
               ylabel="Y Pos", 
               xlim=(0, self._sim_size[0]),
               ylim=(0, self._sim_size[1]),
               aspect=1)
        
        if show: plt.show()