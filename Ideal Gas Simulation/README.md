Web Page Summary: https://isaac-woodard.github.io/projects/idealgassim.html

### Overview
The program is a Python GUI application which simulates the collisions of idealized spherical particles in a three dimensional container. The program was written for use in a student lab in PHYS 309 at Boise State University. 

Features of the program include:
- Modifying the number of particles from 1-500
- Adjusting the temperature of the system by selecting a target temperature
- Changing the simulation speed
- Counting the number of particles within a certain range of speeds
- Displaying an ideal Maxwell-Boltzmann probability distribution for the current temperature

### Dependencies
- Python 3.10+
- Numpy
- matplotlib
- PyQt6
- Numba...
- pyinstaller...

# TODO
Want to refactor the program:
- Use PyQt instead of tkinter
- Create separate classes for the gui and the simulation
- Try to speed up execution with numba
- Create a script to build a standalone executeable with pyinstaller
- Add unit tests