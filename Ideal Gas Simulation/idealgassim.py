import random

import numpy as np
import matplotlib.pyplot as plt
from numba import jit

class Particle:
    def __init__(self, x:float, y:float, z:float, dx:float, dy:float, dz:float, m:int, r:int) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.m = m
        self.r = r
        self.canvasID = None
        self.tagged = None
        
    @property
    def speed(self):
        return (self.dx**2+self.dy**2+self.dz**2)**0.5
    
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
        
        self._sim_depth = radius*6 # Depth in pixels of simulation.
        self._vmax = vmax
        self._vz_max = vmax # Maximum starting depth velocity.
        self._m_per_pix = m_per_pix
        self._vbin_width = vbin_width
        self._vbin_range = vbin_range
        
        self.particles = [] #TODO: Add particles in constructor.
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
                    if x >= p.x-2*p.r and x <= p.x+2*p.r and y >= p.y-2*p.r and y <= p.y+2*p.r and z >= p.z-2*p.r and z <= p.z+2*p.r:
                        overlap = True
                        break
                else:
                    overlap = False 
            p = Particle(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, m=self._mass, r=self._radius)
            nx = int(p.x / self._grid_size)
            ny = int(p.y / self._grid_size)
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

        Args:
            n (int, optional): number of time steps. Defaults to 1.
        """
        if unpaused is True:
            #clears the grid
            global grid
            grid = [[[] for _ in range(height)] for _ in range(length)]
            #update each particle position and the grid
            dt = simSpeed.get() / 50 #parameter to adjust simulation speed
            for p in particleList:
                p.x = p.x + p.dx * dt
                p.y = p.y + p.dy * dt
                p.z = p.z + p.dz * dt #3D Edit Marker
                gridx = int(p.x / gridLength)
                gridy = int(p.y / gridHeight)
                grid[gridx][gridy].append(p)
            #check for collisions with wall and particles and update trajectories if necessary
            for p in particleList:
                if checkWallCollision(p):
                    wallCollision(p)
                else:
                    p2 = checkParticleCollision(p) #store the colliding particle
                    if p2 != p:
                        particleCollision(p, p2)
            drawParticles()
        simulationCanvas.update()
        
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
                speed = mPerPixel * (p.dx**2 + p.dy**2 +p.dz**2)**0.5 #normalizing factor: 200 meters per pixel #3D Edit Marker
                if speed >= low and speed <= high:
                    count += 1
            binCounts[binNum] = (binCounts[binNum]*binTimes+count)/(binTimes+1)
            binNum += 1
        binTimes += 1
        
    #TODO: Refactor
    def check_wall_collision(self, p:Particle) -> bool:
        """Checks the particle for wall collisions.

        Args:
            p (Particle): The particle to check.

        Returns:
            Bool: Whether there is a wall collision.
        """
        if p.x-p.r < 1 or p.x+p.r > simSize or p.y-p.r < 1 or p.y+p.r > simSize or p.z-p.r < 1 or p.z+p.r > simDepth: #3D Edit Marker
            return True
        else:
            return False
        
    #TODO: Refactor
    def wall_collision(self, p:Particle) -> None:
        """Updates particle velocity in case of a wall collision.

        Args:
            p (Particle): The particle to update.
        """
        #bounce the particle
        if p.x-p.r < 1 or p.x+p.r > simSize: #vertical walls
            p.dx = -p.dx
            #make sure the particle isn't pushed behind the wall
            if p.x-p.r < 1:
                p.x = p.r + 1
            if p.x+p.r > simSize:
                p.x = simSize - 2 - p.r
        if p.y-p.r < 1 or p.y+p.r > simSize: #horizontal walls
            p.dy = -p.dy
            #make sure the particle isn't pushed behind the wall
            if p.y-p.r < 1:
                p.y = p.r + 1
            if p.y+p.r > simSize:
                p.y = simSize - 2 - p.r
        #3D Edit Marker
        if p.z-p.r < 1 or p.z+p.r > simDepth: #depth walls
            p.dz = -p.dz
            #make sure the particle isn't pushed behind the wall
            if p.z-p.r < 1:
                p.z = p.r + 1
            if p.z+p.r > simDepth:
                p.z = simDepth - 2 - p.r
            
        #heat or cool the gas based on the target temperature
        if temperature < targetTemp.get():  #heat the gas on collisions with the walls if the temperature is less than the target temperature
            p.dx = p.dx + (targetTemp.get() - temperature)/1000 * p.dx
            p.dy = p.dy + (targetTemp.get() - temperature)/1000 * p.dy
            p.dz = p.dz + (targetTemp.get() - temperature)/1000 * p.dz #3D Edit Marker
            
        if temperature > targetTemp.get(): #cool the gas on collisions with the walls if the temperature is greater than the target temperature
            p.dx = p.dx + (targetTemp.get() - temperature)/1000 * p.dx
            p.dy = p.dy + (targetTemp.get() - temperature)/1000 * p.dy
            p.dz = p.dz + (targetTemp.get() - temperature)/1000 * p.dz #3D Edit Marker
        
    #TODO: Refactor
    def check_particle_collision(self, p:Particle) -> Particle:
        """Checks for particle collisions. Search is confined to neighboring partitions.

        Args:
            p (Particle): Particle to check.

        Returns:
            Particle: The colliding particle if found, otherwise None.
        """
        #find grid position of p
        gridx = int(p.x / gridLength)
        gridy = int(p.y / gridHeight)
        #step through grid spaces surrounding p, as well as the grid space p is in
        for i in range(-1,2):
            for j in range(-1,2):
                if gridx+i < 0 or gridy+j < 0 or gridx+i >= length or gridy+j >= height: #don't check grid spaces beyond the simulation edges (because they don't exist)
                    continue
                #check molecules within the current grid space
                for p2 in grid[gridx+i][gridy+j]:
                #for p2 in particleList:
                    if p == p2:
                        continue
                    mindist = p.r+p2.r #minimum distance between particle centers
                    dx = p.x-p2.x
                    dy = p.y-p2.y
                    dz = p.z-p2.z #3D Edit Marker
                    if (abs(dx) > mindist or abs(dy) > mindist or abs(dz) > mindist): #3D Edit Marker
                        continue
                    dist = (dx**2+dy**2+dz**2)**0.5; #3D Edit Marker
                    if (dist > mindist):
                        continue
                    return p2
        return None
        
    #TODO: Refactor
    def particle_collision(self, p1:Particle, p2:Particle) -> None:
        """Update velocities of colliding particles p1 and p2.

        Args:
            p1 (Particle): The first particle.
            p2 (Particle): The second particle.
        """
        #Algorithm modified from Gas.java (C) 2001 by Paul Falstad, www.falstad.com
     
        #calculate time since the collision theoretically transpired: [(x1-x2)+t(dx1-dx2)]^2 + [(y1-y2)+...]^2 = mindist^2
        dxdif = p1.dx - p2.dx
        xdif  = p1.x - p2.x
        dydif = p1.dy - p2.dy
        ydif  = p1.y - p2.y
        dzdif = p1.dz - p2.dz #3D Edit Marker
        zdif  = p1.z - p2.z #3D Edit Marker
        mindist = p1.r + p2.r #minimum distance
        #a, b and c constants
        a = dxdif**2 + dydif**2 + dzdif**2 #3D Edit Marker
        b = 2*(xdif*dxdif + ydif*dydif + zdif*dzdif) #3D Edit Marker
        c = xdif**2 + ydif**2 + zdif**2 - mindist**2  #3D Edit Marker
        #roots from quadratic formula
        t = (-b - (b**2 - 4*a*c)**0.5)/(2*a)
        t2 = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
        if (abs(t) > abs(t2)):
            t = t2
        
        #move p1 to place and time of collision
        p1.x += t*p1.dx
        p1.y += t*p1.dy
        p1.z += t*p1.dz #3D Edit Marker

        #find the unit vector components between the two particles' centers to treat the collision as one-dimensional
        vx = p1.x - p2.x;
        vy = p1.y - p2.y;
        vz = p1.z - p2.z; #3D Edit Marker
        vxyznorm = (vx**2+vy**2+vz**2)**0.5; #3D Edit Marker
        vxn = vx/vxyznorm;
        vyn = vy/vxyznorm;
        vzn = vz/vxyznorm; #3D Edit Marker

        #treat the particles as one particle and calculate the center of mass velocity
        totmass = p1.m + p2.m
        comdx = (p1.m*p1.dx + p2.m*p2.dx)/totmass
        comdy = (p1.m*p1.dy + p2.m*p2.dy)/totmass
        comdz = (p1.m*p1.dz + p2.m*p2.dz)/totmass #3D Edit Marker

        #find momentum transferred to p1
        pn = (p1.dx - comdx)*vxn + (p1.dy - comdy)*vyn + (p1.dz - comdz)*vzn #3D Edit Marker
        px = 2*vxn*pn
        py = 2*vyn*pn
        pz = 2*vzn*pn #3D Edit Marker

        #transfer momentum to p1
        p1.dx -= px
        p1.dy -= py
        p1.dz -= pz #3D Edit Marker

        #remove momentum from p2
        massratio = p1.m/p2.m
        p2.dx += px*massratio
        p2.dy += py*massratio
        p2.dz += pz*massratio #3D Edit Marker
        
        #End Falstad algorithm
        #go to http://exploratoria.github.io/exhibits/mechanics/elastic-collisions-in-3d/ for more details on the theory
        
    #TODO: Refactor
    def thermometer(self) -> float:
        """Calculate the temperature of the gas based on the Maxwell Boltzmann distribution.

        Returns:
            float: The temperature in Kelvin.
        """
        totalKE = 0
        for p in particleList:
            speed = mPerPixel * p.getSpeed() #200 meters per pixel
            totalKE += 0.5*p.m*1.66e-27*speed**2 #proton mass 1.67e-27 kg
        aveKE = totalKE / numParticles.get() #average the kinetic energy
        global temperature
        temperature = aveKE / 1.38e-23 / 1.5 #find the temperature
        #update temperature value for label
        tempString.set(f"T = {temperature:.1f}")