/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** multigrid.h **************************

#include <stdio.h>            
#include <stdlib.h>           
#include <math.h>             
#include <string.h>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace std;

#include <globals.h>
#include <fileio.h>
#include <hdf5write.h>
#include <array3d.h>

class MultiGrid
{
 public:

  double w;		// Successive Over-Relaxation factor
  int ncycle;		// Number of SOR sysles at each resolution
  int iterations;		// Number of VCycles

  int ScaleFactor;       // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 1.0 micron, 96 grids in the z-direction
  double PixelSize;      // Pixel size in microns
  int GridsPerPixel;     // Number of grids per pixel at ScaleFactor = 1
  int Nx;                // Number of grids in x at ScaleFactor = 1
  int Ny;                // Number of grids in y at ScaleFactor = 1
  int Nz;                // Number of grids in z at ScaleFactor = 1
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
  int nsteps;            // Number of steps in Multigrid

  // Voltages and Charges
  double Vbb;		// Back bias
  double Vparallel_lo;	// Parallel Low Voltage
  double Vparallel_hi;	// Parallel High Voltage
  double Vserial_lo;	// Serial Low Voltage
  double Vaverage;       // Average volatge on bottom
  double Vchannelstop;   //  Voltage of undepleted channel stop
  double Vscupper;       // Scupper voltage

  double BackgroundDoping; 	// Background doping
  double ChannelStopDoping;	// Channel Stop doping
  double ChannelStopDepth;     	// Channel stop depth in microns
  double ChannelStopWidth;     	// Channel stop width in microns
  double ChannelDoping;		// Channel doping
  double ChannelDepth;		// Channel depth in microns
  int UndepletedChannelStop;	// 0 = No undepleted Region, 1 = undepleted Region

// Pixel Regions

  int NumberofPixelRegions;	  	  // 1
  double** PixelRegionLowerLeft;	  // 
  double** PixelRegionUpperRight;	  //
  int* NumberofFilledWells;		  //
  int** CollectedCharge;		  // Collected charge in e-
  double*** FilledPixelCoords;            // (x,y) coords of pixel center
  int DistributedCharge;                   // 1=Charge distributed over well, 0=Point charge in center
  
// Added dipoles

  int NumberofDipoles;                  // Number of dipoles
  double** DipoleCoords;                // X,Y coords of dipole center
  double* DipoleZLocation;               // Z location of dipole charge
  int* DipoleCharge;                     // Number of electrons in dipole
  int NumberofDipoleImages;             // Number of image pairs to calculate
  
// Constant Voltage Regions

  int NumberofFixedRegions;
  double** FixedRegionLowerLeft;
  double** FixedRegionUpperRight;
  double* FixedRegionVoltage;
  double* FixedRegionDoping;
  double* FixedRegionChargeDepth;

  // Pixel Boundary Tests

  int LogEField;
  int LogPixels;
  int LogPixelPaths;
  double* PixelBoundaryLowerLeft;
  double* PixelBoundaryUpperRight;
  double* PixelBoundaryStepSize;

  string outputfilebase; // Output filename base
  string outputfiledir; // Output filename directory

  Array3D** phi;      // Phi arrays
  Array3D** rho;      // Rho arrays
  Array3D** E;
  
  MultiGrid() {};
  MultiGrid(string);
  ~MultiGrid();
 
  void ReadConfigurationFile(string);
  void BuildArrays(Array3D**);
  void SetInitialConditions(Array3D*, Array3D*);
  void SOR(Array3D*, Array3D*, double);
  void RedBlack(Array3D*, Array3D*);
  double Error(Array3D*, Array3D*);
  void Restrict(Array3D*, Array3D*, Array3D*, Array3D*);
  void Prolongate(Array3D*, Array3D*);
  void VCycle(Array3D**, Array3D**, double, int, int);
  void WriteOutputFile(string, string, string, Array3D*);
  void Gradient(Array3D*, Array3D**);
  void Trace(Array3D*, Array3D**, double*);
  void TracePaths(Array3D*, Array3D**, double*, ofstream&);
  void TraceGrid(Array3D*, Array3D**);
  void AddDipolePotentials(Array3D*);
};

