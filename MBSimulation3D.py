'''
Created on Dec 14, 2018
Last Updated: Aug 26, 2019
@author: Isaac Woodard
'''

from tkinter import *
from tkinter import ttk
import random
import matplotlib.pyplot as plt
import numpy as np

#Notes: There is a bug that crashes the program if the particles move too fast. Currently, this is avoided by limiting the temperature to 500 K.
#Notes: The simulation seems to run faster for lower particle counts.


#########################
#Adjustable Parameters
#########################
simSize = 700            #width and height in pixels of the simulation frame
defaultNumPars = 250     #starting number of particles
defMass = 4              #default mass of particles in units of amu
defRadius = 5            #default radius of particles in pixels

####################
#Other Parameters
####################
simDepth = defRadius*6      #depth in pixels of simulation
maxVelocity = 6.8           #maximum starting horizontal and vertical velocity of all particles in terms of pixels per timestep
maxzVelocity = maxVelocity  #maximum starting depth velocity of all particles in terms of pixels per timestep
mPerPixel = 200             #meters per pixel
#the following parameters are used for aggregated data collection for debugging
binRange = 200              #"width" of the velocity bins
maxRange = 3000             #max range of velocity binning

########################################################
#Support Functions for GUI (Graphical User Interface)
########################################################

#initialize variables
lowspeed = 0
highspeed = 0
unpaused = True
temperature = 0
binCounts = [0 for _ in range(int(maxRange/binRange))]
binTimes = 0

#used to pause and unpause the program when the "pause" button is clicked
def pause():
    global unpaused
    if unpaused is True:
        unpaused = False
    else:
        unpaused = True
    
#used to reset the simulation by removing all particles and creating new ones
def reset():
    removeParticles(len(particleList))
    addParticles(defaultNumPars)
    numParticles.set(defaultNumPars)
    numParticlesString.set("N = 500")
    targetTempString.set("T = 294")
    targetTemp.set(294)
    simSpeed.set(50)
    global binCounts, binTimes
    binCounts = [0 for _ in range(int(maxRange/binRange))]
    binTimes = 0
    
    
#used to reset the aggregated data for debugging
def resetBins():
    global binCounts, binTimes
    binCounts = [0 for _ in range(int(maxRange/binRange))]
    binTimes = 0

#used to adjust the number of particles. Also updates the label for the number of particles.
def adjustParticleCount(*args):
    if len(particleList) < numParticles.get():
        a = numParticles.get() - len(particleList)
        addParticles(a)
    if len(particleList) > numParticles.get():
        b = len(particleList) - numParticles.get()
        removeParticles(b)
    numParticlesString.set("N =" + str(numParticles.get()))
    
def adjustTargetTemp(*args):
    targetTempString.set("Target T = " + str(targetTemp.get()))
    
def drawParticles():
    #TODO: update particle colors based on speed
    simulationCanvas.delete("all")
    for p in particleList:
        mRad = int(p.r * (0.33 + 0.67*p.z/simDepth)) #modifies the size of the circle drawn based on the particle's depth
        if p.tagged is 0:
            p.canvasID = simulationCanvas.create_oval(p.x-mRad, p.y-mRad, p.x+mRad, p.y+mRad, fill="blue");
        if p.tagged is 1:
            p.canvasID = simulationCanvas.create_oval(p.x-mRad, p.y-mRad, p.x+mRad, p.y+mRad, fill="green");
    
#used to count the number of particles in a certain speed range
def countSpeeds():   
    try:
        count = 0
        #untag all particles
        for p in particleList:
            p.tagged = 0
        #get high and low speed values from text fields
        low = int(lowspeed.get())
        high = int(highspeed.get())
        #step through particles and find speeds
        for p in particleList:
            speed = mPerPixel * (p.dx**2 + p.dy**2 +p.dz**2)**0.5 #normalizing factor: 200 meters per pixel #3D Edit Marker
            if speed >= low and speed <= high:
                count += 1
                p.tagged = 1
        #update label
        countString.set(str(count))
        if unpaused is False:
            drawParticles()
    except ValueError:
        print("Speed range must be given with integer values.")
        
#used to get data for the theoretical Maxwell-Boltzmann Distribution Plot
def mbPlot():
    T = temperature #temperature in Kelvin
    m = defMass*1.66e-27 #mass in kilograms
    k = 1.381e-23 #Boltzmann's constant in Joules per Kelvin
    
    # Data for plotting theoretical Maxwell-Boltzmann distribution for two dimensions
    v = np.arange(0.0, maxRange, 0.1) #x values, velocity
    dNdv = numParticles.get()*(m/2/np.pi/k/T)**1.5*4*np.pi*v**2*np.exp(-m*v**2/2/k/T) #y values, probability from Boltzmann distribution times number of particles #3D Edit Marker
    plt.plot(v, dNdv)
    plt.xlabel('Speed (m/s)')
    plt.ylabel('dN/dv')
    plt.title(f"MB Distribution: {temperature:.1f} K")
    plt.show(block=False)
    
#used to print the average bin count for each speed range
def binCount():
    print("Samples: ",binTimes)
    j=0
    particleSum = 0
    for i in binCounts:
        print(binRange/2+j*binRange,": ",'%.3f'%(i))
        particleSum += i
        j+=1
    print("Particles: ",particleSum)
    print("---------------\n")
    
###################
#Code for GUI
###################

root = Tk()
root.title("Particle Simulation")
root.resizable(False, False)

mainFrame = ttk.Frame(root, padding="2 1 2 1")
mainFrame.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=2)
root.rowconfigure(0, weight=1)

simulationCanvas = Canvas(mainFrame, width=simSize, height=simSize, background="black", borderwidth=2, relief="sunken")
simulationCanvas.grid(column=0, row=0, sticky=(N, W, E, S))

interfaceFrame = ttk.Frame(mainFrame, padding="10 5 10 1")
interfaceFrame.grid(column=1, row=0, sticky=(N, W, E, S))

#TODO: make a combobox with different options for initial parameters such as particle size, mass and initial speeds
#parameters = StringVar()
#parameterBox = ttk.Combobox(interfaceFrame, textvariable=parameters, state="readonly").grid(column=0, row=0)
#parameterBox['parameters'] = ('Standard_Settings', 'B')
#parameterBox.current(0)

controlLabel = ttk.Label(interfaceFrame, text="Controls", font="bold").grid(column=0, row=1, pady="20 5")
controlFrame = ttk.Frame(interfaceFrame, borderwidth=3, relief="groove")

prFrame = ttk.Frame(controlFrame, padding="5 5 5 5")
pauseButton = ttk.Button(prFrame, text="Pause", command=pause).grid(column=0, row=0, padx=5)
resetButton = ttk.Button(prFrame, text="Reset", command=reset).grid(column=1, row=0, padx=5)
prFrame.grid(column=0, row=0)

numParticles = IntVar()
particleLabel = ttk.Label(controlFrame, text="Particle Count").grid(column=0, row=1, pady=[10,0])
particleSlider = ttk.Scale(controlFrame, orient=HORIZONTAL, length=200, from_=1.0, to=500.0, variable=numParticles, command=adjustParticleCount).grid(column=0, row=2, padx=30)
numParticles.set(defaultNumPars)

targetTemp = IntVar()
reservoirLabel = ttk.Label(controlFrame, text="Target Temperature").grid(column=0, row=3, pady=[10,0])
reservoirSlider = ttk.Scale(controlFrame, orient=HORIZONTAL, length=200, from_=0, to=500, variable=targetTemp, command=adjustTargetTemp).grid(column=0, row=4)
targetTemp.set(294)

simSpeed = IntVar()
speedLabel = ttk.Label(controlFrame, text="Simulation Speed").grid(column=0, row=5, pady=[10,0])
speedSlider = ttk.Scale(controlFrame, orient=HORIZONTAL, length=200, from_=0.0, to=100, variable=simSpeed).grid(column=0, row=6)
simSpeed.set(50) 

controlFrame.grid(column=0, row=2)

statsLabel = ttk.Label(interfaceFrame, text="Data", font="bold").grid(column=0, row=3, pady="20 5")
statsFrame = ttk.Frame(interfaceFrame, borderwidth=3, relief="groove")

varFrame = ttk.Frame(statsFrame, padding="5 5 5 5")
numParticlesString = StringVar()
numParticlesString.set("N = 500")
numberLabel = ttk.Label(varFrame, textvariable=numParticlesString).grid(column=0, row=0, padx=5)
tempString = StringVar()
tempString.set("T = ---")
tempLabel = ttk.Label(varFrame, textvariable=tempString).grid(column=1, row=0, padx=5)
targetTempString = StringVar()
targetTempString.set("Target T = 294")
targetTempLabel = ttk.Label(varFrame, textvariable=targetTempString).grid(column=2, row=0, padx=5)
varFrame.grid(column=0, row=0)

searchFrame = ttk.Frame(statsFrame, padding="5 5 5 5")
lowspeed = StringVar()
lowLabel = ttk.Label(searchFrame, text="Low Speed (m/s)").grid(column=0, row=0)
lowEntry = ttk.Entry(searchFrame, textvariable=lowspeed).grid(column=0, row=1)
highspeed = StringVar()
highLabel = ttk.Label(searchFrame, text="High Speed (m/s)").grid(column=1, row=0)
highEntry = ttk.Entry(searchFrame, textvariable=highspeed).grid(column=1, row=1)
countButton = ttk.Button(searchFrame, text="Get Bin Count", command=countSpeeds).grid(column=0, row=2)
countString = StringVar()
countString.set("---")
countLabel = ttk.Label(searchFrame, textvariable=countString).grid(column=1, row=2)
searchFrame.grid(column=0, row=1)

graphButton = ttk.Button(statsFrame, text="Get Ideal MB Graph", command=mbPlot).grid(column=0, row=3, pady=5)
#The following buttons are used to output large amounts of data for bug testing
#binButton = ttk.Button(statsFrame, text="Get Average Bin Count", command=binCount).grid(column=0, row=4, pady=5)
#binresetButton = ttk.Button(statsFrame, text="Reset Bins", command=resetBins).grid(column=0, row=5, pady=5)

statsFrame.grid(column=0, row=4)

#spin = Spinbox(parent, from_=1.0, to=100.0, textvariable=spinval)

################################
#Particle Simulation
################################

particleList = [] #the complete list of particles
gridLength = 50 #the pixel length of each grid spacing
gridHeight = 50 #the pixel height of each grid spacing
length = int(simSize / gridLength + 1) #number of horizontal grid spacings
height = int(simSize / gridHeight + 1) #number of vertical grid spacings
grid = [[[] for _ in range(height)] for _ in range(length)] #the grid partition of particles

class Particle:
    x = -1 #horizontal position of particle center
    y = -1 #vertical position of particle center
    z = -1 #depth position of particle center #3D Edit Marker
    dx = 0.0 #horizontal velocity, measured in pixels
    dy = 0.0 #vertical velocity, measured in pixels
    dz = 0.0 #depth velocity, measured in pixels #3D Edit Marker
    m = 0 #mass
    r = 0 #radius
    canvasID = -1 #canvasID number of canvas object
    tagged = 0 #1 if the particle was one of the last counted by the user
    def getSpeed(self):
        return (self.dx**2+self.dy**2+self.dz**2)**0.5 #3D Edit Marker
    def getKE(self):
        s = self.getSpeed
        return 0.5*self.m*s**2
    def getMomentum(self):
        s = self.getSpeed
        return self.m*s

#used both to begin the simulation and add particles in the middle of it
def addParticles(n):
    for i in range(n):
        global IDcounter
        #create particle
        p = Particle()
        p.m = defMass
        p.r = defRadius
        p.dx = random.uniform(-maxVelocity, maxVelocity)
        p.dy = random.uniform(-maxVelocity, maxVelocity)
        p.dz = random.uniform(-maxzVelocity, maxzVelocity) #3D Edit Marker
        #place particle & check against other positions
        if len(particleList) is 0:
            x = random.randint(1+p.r, simSize-p.r)
            y = random.randint(1+p.r, simSize-p.r)
            z = random.randint(1+p.r, simDepth-p.r) #3D Edit Marker
        else: 
            check = True
            while check: #while there might be overlap between particles
                x = random.randint(1+p.r, simSize-p.r)
                y = random.randint(1+p.r, simSize-p.r)
                z = random.randint(1+p.r, simDepth-p.r) #3D Edit Marker
                for i in particleList: #check overlapping with each particle
                    if x >= i.x-2*i.r and x <= i.x+2*i.r and y >= i.y-2*i.r and y <= i.y+2*i.r and z >= i.z-2*i.r and z <= i.z+2*i.r: #check horizontal, vertical ranges and depth ranges #3D Edit Marker
                        check = True
                        break
                else:
                    check = False 
        
        #draw particle
        p.x = x
        p.y = y
        p.z = z #3D Edit Marker
        mRad = int(p.r * (0.33 + 0.67*p.z/simDepth)) #modifies the size of the circle drawn based on the particle's depth
        p.canvasID = simulationCanvas.create_oval(p.x-mRad, p.y-mRad, p.x+mRad, p.y+mRad, fill="blue");
        #add particle to grid and master list
        gridx = int(p.x / gridLength)
        gridy = int(p.y / gridHeight)
        grid[gridx][gridy].append(p)
        particleList.append(p) #master list

#used when taking particles out of the system (and maybe to reset)
def removeParticles(n):
    for _ in range(n):
        p = particleList.pop()
        simulationCanvas.delete(p.canvasID)

#used to advance the simulation one time-step, redrawing ovals and checking for collisions as necessary
def timeStep():  
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
        
#used to check each particle's movement for wall collisions on the following time-step
def checkWallCollision(p):
    if p.x-p.r < 1 or p.x+p.r > simSize or p.y-p.r < 1 or p.y+p.r > simSize or p.z-p.r < 1 or p.z+p.r > simDepth: #3D Edit Marker
        return True
    else:
        return False

#used to update particle velocity in case of a wall collision
def wallCollision(p):
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

#used to check for particle collisions using the grid system 
def checkParticleCollision(p):
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
                return p2 #if a collision is found
    return p #if no collision is found

#used to update velocities in case of two particles colliding
def particleCollision(p1, p2):
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
    
#used to calculate the temperature of the particles based on the Maxwell Boltzmann distribution and update the temperature label
def thermometer():
    #find total kinetic energy
    totalKE = 0
    for p in particleList:
        speed = mPerPixel * p.getSpeed() #200 meters per pixel
        totalKE += 0.5*p.m*1.66e-27*speed**2 #proton mass 1.67e-27 kg
    aveKE = totalKE / numParticles.get() #average the kinetic energy
    global temperature
    temperature = aveKE / 1.38e-23 / 1.5 #find the temperature
    #update temperature value for label
    tempString.set(f"T = {temperature:.1f}")
    
#######################
#Code to Run Program
#######################
    
addParticles(numParticles.get())
while True: #loops the simulation
    try:
        root.after(1, timeStep())
    except: #TODO: pinpoint '_tkinter.TclError: invalid command name ".!frame.!canvas"' and/or find a solution that doesn't involve ignoring the error
        quit()
    thermometer()
root.mainloop()