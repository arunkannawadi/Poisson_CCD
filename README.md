# Poisson_CCD
Poisson solver for LSST CCDs
Description of stand-alone Poisson solver.
Craig Lage - UC Davis - 18-Feb-15

Description: This code is a simple grid-based Poisson's equation solver intended to simulate pixel distortion effects in thick fully-depleted CCD's to be used in the LSST digital camera.  The code builds a 3D rectilinear grid to represent a portion of the CCD, assigns the appropriate charge densities and applied potentials, then solves Poisson's equation using multi-grid methods.  A 256^3 grid, which is adequate for most purposes solves in less than one minute on a typical laptop.  Plotting routines are available to plot the potentials, E-Fields, Pixel shapes, and electron paths.  A brief description of the code and some samples of the output are in the file docs/Poisson.pdf.  Below is a basic description of how to install and run the code:

Installing: Read the INSTALL file.

Running:  The basic syntax is:

./Poisson <configuration file>

If you have compiled the code as described in INSTALL, the executable will be in src/Poisson.  There are four example configuration files, pixels_256.cfg, edge_256.cfg, pixels_512.cfg, and edge_512.cfg, in the examples folder.  For example, to run the pixels_256.cfg file, cd into the Poisson folder and type:

src/Poisson examples/pixels_256.cfg

This will generate a series of files in the data folder (or wherever you tell it to put the files in the .cfg file).
Most of the settings in the .cfg file should be relatively self-explanatory.

Plotting:  A Python plotting routine is available that makes a series of plots.  You will need to install the h5py module so that Python can read the hdf5 files.  Run this program by typing:

python Poisson_Plots.py <configuration file>

The plotting routine should run OK with the four example .cfg files, but if you do something else, you will probably need to customize the Poisson_Plots.py code.

Hopefully you find this useful.  Comments and questions are encouraged and should be addressed to: cslage@ucdavis.edu

