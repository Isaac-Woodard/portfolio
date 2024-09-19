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
        ...
        
    def reset():
        ...
        
    def reset_bins():
        ...
        
    def update_nparticles():
        ...
        
    def update_target_temp():
        ...
        
    def draw_particles():
        ...
        
    def count_speeds():
        ...
        
    def plot_MB():
        ...
        
    def bin_count():
        ...
        
if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)
    app = IdealGasGUI()
    sys.exit(qapp.exec())