Web Page Summary: https://isaac-woodard.github.io/projects/dlasim.html

# Overview
Diffusion-Limited Aggregation (DLA) is the growth of a cluster, such as a crystal, one particle at a time. The particles move by diffusion through the surrounding medium before making contact with the seed for the structure. The resulting cluster is fractal in nature.

This program creates 2D DLA clusters to model 4 different physical phenomena. The simulation is written in Python and organized in a Jupyter Notebook. The first 8 cells create the clusters for each phenomenon. The second to last cell defines the function used to create the clusters. 

The final cell contains a function for calculating the box-counting dimension (fractal dimension) of a cluster (found at https://github.com/rougier/numpy-100).

# Usage
Go to https://jupyter.org/ and either download and install Jupyter Notebook for desktop or select the "Try it in your browser" button. If using the browser, create a Classic Notebook and select File > Open. Upload the notebook file to the binder.

To create a new cluster, either copy the contents of the first two cells into two new cells or modify them directly. The first three variables define the parameters for the simulation. Press ctrl+enter or click the play button to evaluate a cell.
