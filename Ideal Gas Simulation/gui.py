import sys

from PyQt6 import QtWidgets, uic, QtCore

import idealgassim as gas

class IdealGasGUI(QtWidgets.QMainWindow):
    """Defines a GUI for the ideal gas simulation."""
    
    def __init__(self, show:bool=True) -> None:
        """Initializes the GUI.

        Args:
            show (bool, optional): Whether to show the UI. Used for testing.
        """
        super(IdealGasGUI, self).__init__()
        uic.loadUi("GUI.ui", self)
        
        #TODO: Assign to self
        lowspeed = 0
        highspeed = 0
        unpaused = True
        temperature = 0
        binCounts = [0 for _ in range(int(maxRange/binRange))]
        binTimes = 0
        
        #TODO: Add Listeners
        
        if show: self.show()
        
    def pause():
        #used to pause and unpause the program when the "pause" button is clicked
        global unpaused
        if unpaused is True:
            unpaused = False
        else:
            unpaused = True
        
    def reset():
        #used to reset the simulation by removing all particles and creating new ones
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
        
    def reset_bins():
        #used to reset the aggregated data for debugging
        global binCounts, binTimes
        binCounts = [0 for _ in range(int(maxRange/binRange))]
        binTimes = 0
        
    def update_nparticles():
        #used to adjust the number of particles. Also updates the label for the number of particles.
        if len(particleList) < numParticles.get():
            a = numParticles.get() - len(particleList)
            addParticles(a)
        if len(particleList) > numParticles.get():
            b = len(particleList) - numParticles.get()
            removeParticles(b)
        numParticlesString.set("N =" + str(numParticles.get()))
        
    def update_target_temp():
        targetTempString.set("Target T = " + str(targetTemp.get()))
        
    def draw_particles():
        simulationCanvas.delete("all")
        for p in particleList:
            mRad = int(p.r * (0.33 + 0.67*p.z/simDepth)) #modifies the size of the circle drawn based on the particle's depth
            if p.tagged is 0:
                p.canvasID = simulationCanvas.create_oval(p.x-mRad, p.y-mRad, p.x+mRad, p.y+mRad, fill="blue");
            if p.tagged is 1:
                p.canvasID = simulationCanvas.create_oval(p.x-mRad, p.y-mRad, p.x+mRad, p.y+mRad, fill="green");
       
    #used to count the number of particles in a certain speed range 
    def count_speeds():
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
        
    def plot_MB():
        #used to get data for the theoretical Maxwell-Boltzmann Distribution Plot
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
        
    def bin_count():
        #used to print the average bin count for each speed range
        print("Samples: ",binTimes)
        j=0
        particleSum = 0
        for i in binCounts:
            print(binRange/2+j*binRange,": ",'%.3f'%(i))
            particleSum += i
            j+=1
        print("Particles: ",particleSum)
        print("---------------\n")
        
if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)
    app = IdealGasGUI()
    sys.exit(qapp.exec())