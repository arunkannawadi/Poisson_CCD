/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** multigrid.cpp **************************

#include "multigrid.h"

MultiGrid::MultiGrid(string inname) //Constructor                                                                                            
{
  // This reads in the data from the poisson.cfg
  // file, sets the initial conditions, and solves
  // Poisson's equation using SOR and Multi-Grid methods

  // First we read in the configuration information
  ReadConfigurationFile(inname);

  // Then, we build the multigrid arrays and set the initial conditions
  nsteps = 5 + (int)(log2(ScaleFactor));
  // nsteps is the number of reduction steps in the VCycle.
  // This gives 4 grid cells in the z-direction at the coarsest scale
  phi = new Array3D*[nsteps+1];
  rho = new Array3D*[nsteps+1];
  BuildArrays(phi);
  BuildArrays(rho);
  SetInitialConditions(phi[0], rho[0]);
  printf("Finished Setting ICs\n");
  fflush(stdout);
  // Now we cycle through the VCycles to solve Poisson's equation
  int n;
  for (n=0; n<iterations; n++)
    {
      VCycle(phi, rho, w, nsteps, ncycle);
    }

  // Add dipoles if there are any:
  if (NumberofDipoles>0)
    {
      AddDipolePotentials(phi[0]);
    }
  printf("Finished adding dipole potentials.\n");
  fflush(stdout);      
  // Now, we write out the potential and charge density
  WriteOutputFile(outputfiledir, outputfilebase, "phi", phi[0]);
  WriteOutputFile(outputfiledir, outputfilebase, "rho", rho[0]);

  // Next we calculate the E fields (if asked)
  if (LogEField == 1)
    {
      E = new Array3D*[3];
      Gradient(phi[0], E);
      WriteOutputFile(outputfiledir, outputfilebase, "Ex", E[0]);
      WriteOutputFile(outputfiledir, outputfilebase, "Ey", E[1]);
      WriteOutputFile(outputfiledir, outputfilebase, "Ez", E[2]);
    }
  // Last we calculate the pixel boundaries (if asked, and if E-Feilds exist)
  if (LogPixels == 1  && LogEField == 1)
    {
      TraceGrid(phi[0], E);
    }
  return;
}

MultiGrid::~MultiGrid() //Destructor                                                                                            
{
  int n;
  for (n=0; n<nsteps+1; n++)
    {
      delete phi[n];
      delete rho[n];
    }
  if (LogEField == 1)
    {
      for (n=0; n<3; n++)
	{
	  delete E[n];
	}
      delete[] E;
    }
  delete[] phi;
  delete[] rho;
  return;
}

void MultiGrid::ReadConfigurationFile(string inname)
{
  // Poisson solver constants
  w = GetDoubleParam(inname, "w", 1.9);			// Successive Over-Relaxation factor
  ncycle = GetIntParam(inname, "ncycle", 100);		// Number of SOR sysles at each resolution
  iterations =  GetIntParam(inname, "iterations", 3);	// Number of VCycles

  // Overall setup
  
  ScaleFactor =  GetIntParam(inname, "ScaleFactor", 1);     // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 5/6 micron, 128 grids in the z-direction
  PixelSize = GetDoubleParam(inname, "PixelSize", 10.0);    // Pixel size in microns
  GridsPerPixel = GetIntParam(inname, "GridsPerPixel", 12); // Grids per pixel at ScaleFactor = 1
  GridsPerPixel = GridsPerPixel * ScaleFactor;
  Nx = GetIntParam(inname, "Nx", 128);                // Number of grids in x at ScaleFactor = 1
  Nx = Nx * ScaleFactor;
  Ny = GetIntParam(inname, "Ny", 128);                // Number of grids in y at ScaleFactor = 1
  Ny = Ny * ScaleFactor;
  Nz = 128;                                           // Number of grids in z at ScaleFactor = 1
  Nz = Nz * ScaleFactor;
  
  // Voltages and Charges
  Vbb = GetDoubleParam(inname, "Vbb", -50.0);		        // Back bias
  Vparallel_lo = GetDoubleParam(inname, "Vparallel_lo", -8.0);	// Parallel Low Voltage
  Vparallel_hi = GetDoubleParam(inname, "Vparallel_hi", 4.0);	// Parallel High Voltage
  Vserial_lo = GetDoubleParam(inname, "Vserial_lo", -6.0);	// Serial Low Voltage
  Vaverage = (8.0 * Vparallel_lo + 4.0 * Vparallel_hi) / 12.0;
  Vchannelstop = GetDoubleParam(inname, "Vchannelstop", 0.0);
  Vscupper = GetDoubleParam(inname, "Vscupper", 5.0);

  BackgroundDoping = GetDoubleParam(inname, "BackgroundDoping", -1.0E12);
  ChannelStopDoping = GetDoubleParam(inname, "ChannelStopDoping", -1.0E12);
  ChannelStopDepth = GetDoubleParam(inname, "ChannelStopDepth", 1.0);
  ChannelStopWidth = GetDoubleParam(inname, "ChannelStopWidth", 1.0);
  ChannelDoping = GetDoubleParam(inname, "ChannelDoping", -5.0E11);
  ChannelDepth = GetDoubleParam(inname, "ChannelDepth", 1.0);
  UndepletedChannelStop = GetIntParam(inname, "UndepletedChannelStop", 0);

  // Pixel Regions
  int i, j, k;
  string regionnum, fillednum;
  NumberofPixelRegions = GetIntParam(inname, "NumberofPixelRegions", 0);
  PixelRegionLowerLeft = new double*[NumberofPixelRegions];
  PixelRegionUpperRight = new double*[NumberofPixelRegions];
  NumberofFilledWells = new int[NumberofPixelRegions];
  DistributedCharge = GetIntParam(inname, "DistributedCharge", 1);
  FilledPixelCoords = new double**[NumberofPixelRegions];
  CollectedCharge = new int*[NumberofPixelRegions];

  for (i=0; i<NumberofPixelRegions; i++)
    {
      PixelRegionLowerLeft[i] = new double[2];
      PixelRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  PixelRegionLowerLeft[i][j] = 0.0;
	  PixelRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofPixelRegions; i++)
    {
      regionnum = boost::lexical_cast<std::string>(i);
      PixelRegionLowerLeft[i] = GetDoubleList(inname, "PixelRegionLowerLeft_"+regionnum, 2, PixelRegionLowerLeft[i]);
      PixelRegionUpperRight[i] = GetDoubleList(inname, "PixelRegionUpperRight_"+regionnum, 2, PixelRegionUpperRight[i]);
      NumberofFilledWells[i] = GetIntParam(inname, "NumberofFilledWells_"+regionnum, 0);
      CollectedCharge[i] = new int[NumberofFilledWells[i]];
      FilledPixelCoords[i] = new double*[NumberofFilledWells[i]];
      for (j=0; j<NumberofFilledWells[i]; j++)
	{
	  fillednum = boost::lexical_cast<std::string>(j);
	  CollectedCharge[i][j] = GetIntParam(inname,"CollectedCharge_"+regionnum+"_"+fillednum,0);
	  FilledPixelCoords[i][j] = new double[2];
	  for (k=0; k<2; k++)
	    {
	      FilledPixelCoords[i][j][k] = 0.0;
	    }
	  FilledPixelCoords[i][j] = GetDoubleList(inname, "FilledPixelCoords_"+regionnum+"_"+fillednum, 2, FilledPixelCoords[i][j]);
	}
    }

  // Dipoles added

  string dipolenum;
  NumberofDipoles = GetIntParam(inname, "NumberofDipoles", 0);
  NumberofDipoleImages = GetIntParam(inname, "NumberofDipoleImages", 5);
  DipoleCoords = new double*[NumberofDipoles];
  DipoleZLocation = new double[NumberofDipoles];
  DipoleCharge = new int[NumberofDipoles];

  for (i=0; i<NumberofDipoles; i++)
    {
      DipoleCoords[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  DipoleCoords[i][j] = 0.0;
	}
    }
  for (i=0; i<NumberofDipoles; i++)
    {
      dipolenum = boost::lexical_cast<std::string>(i);
      DipoleCoords[i] = GetDoubleList(inname, "DipoleCoords_"+dipolenum, 2, DipoleCoords[i]);
      DipoleZLocation[i] = GetDoubleParam(inname,"DipoleZLocation_"+dipolenum,0.5);
      DipoleCharge[i] = GetIntParam(inname,"DipoleCharge_"+dipolenum,0);
    }
  
  // Fixed Voltage Regions
  NumberofFixedRegions = GetIntParam(inname, "NumberofFixedRegions", 0);
  FixedRegionLowerLeft = new double*[NumberofFixedRegions];
  FixedRegionUpperRight = new double*[NumberofFixedRegions];
  FixedRegionVoltage = new double[NumberofFixedRegions];
  FixedRegionDoping = new double[NumberofFixedRegions];
  FixedRegionChargeDepth = new double[NumberofFixedRegions];

  for (i=0; i<NumberofFixedRegions; i++)
    {
      FixedRegionLowerLeft[i] = new double[2];
      FixedRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  FixedRegionLowerLeft[i][j] = 0.0;
	  FixedRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofFixedRegions; i++)
    {
      regionnum = boost::lexical_cast<std::string>(i);
      FixedRegionLowerLeft[i] = GetDoubleList(inname, "FixedRegionLowerLeft_"+regionnum, 2, FixedRegionLowerLeft[i]);
      FixedRegionUpperRight[i] = GetDoubleList(inname, "FixedRegionUpperRight_"+regionnum, 2, FixedRegionUpperRight[i]);
      FixedRegionVoltage[i] = GetDoubleParam(inname, "FixedRegionVoltage_"+regionnum,0.0);
      FixedRegionDoping[i] = GetDoubleParam(inname, "FixedRegionDoping_"+regionnum,0.0);
      FixedRegionChargeDepth[i] = GetDoubleParam(inname, "FixedRegionChargeDepth_"+regionnum,0.0);
    }


  // Pixel Boundary Tests

  LogEField = GetIntParam(inname, "LogEField", 0);
  LogPixels = GetIntParam(inname, "LogPixels", 0);  
  PixelBoundaryLowerLeft = new double(2);
  PixelBoundaryUpperRight = new double(2);
  PixelBoundaryStepSize = new double(2);
  for (j=0; j<2; j++)
    {
      PixelBoundaryLowerLeft[j] = 0.0;
      PixelBoundaryUpperRight[j] = 100.0;
      PixelBoundaryStepSize[j] = 1.0;
    }
  PixelBoundaryLowerLeft = GetDoubleList(inname, "PixelBoundaryLowerLeft", 2, PixelBoundaryLowerLeft);
  PixelBoundaryUpperRight = GetDoubleList(inname, "PixelBoundaryUpperRight", 2, PixelBoundaryUpperRight);
  PixelBoundaryStepSize = GetDoubleList(inname, "PixelBoundaryStepSize", 2, PixelBoundaryStepSize);
  LogPixelPaths = GetIntParam(inname, "LogPixelPaths", 0);
  
  outputfilebase  = GetStringParam(inname,"outputfilebase", "Test"); //Output filename base
  outputfiledir  = GetStringParam(inname,"outputfiledir", "data"); //Output filename directory
  return;
}

void MultiGrid::BuildArrays(Array3D** array)
{
  // Builds the multigrid arrays
  int nx, ny, nz, nxx, nyy, nzz, n;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, CenterShift;
  nxx = Nx + 1;
  nyy = Ny + 1;
  nzz = Nz + 1;
  dx = PixelSize / (double)GridsPerPixel;
  dy = dx;
  dz = dx;
  CenterShift = PixelSize / (2.0 * (double)GridsPerPixel);
  Xmin = -dx / 2.0 + CenterShift;
  Ymin = -dy / 2.0 + CenterShift;
  Zmin = -dz / 2.0 + CenterShift;
  Xmax = dx * ((double)nxx - 1.0 / 2.0) + CenterShift;
  Ymax = dy * ((double)nyy - 1.0 / 2.0) + CenterShift;
  Zmax = dz * ((double)nzz - 1.0 / 2.0) + CenterShift;
  array[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz);

  for (n=1; n<nsteps+1; n++)
    {
      nx = (array[0]->nx - 1) / (int)pow(2,n) + 1;
      ny = (array[0]->ny - 1) / (int)pow(2,n) + 1;
      nz = (array[0]->nz - 1) / (int)pow(2,n) + 1;
      dx = array[0]->dx * (int)pow(2,n);
      dy = array[0]->dy * (int)pow(2,n);
      dz = array[0]->dz * (int)pow(2,n);
      xmin = array[0]->xmin + array[0]->dx / 2.0 - dx / 2.0;
      ymin = array[0]->ymin + array[0]->dy / 2.0 - dy / 2.0;
      zmin = array[0]->zmin + array[0]->dz / 2.0 - dz / 2.0; 
      xmax = array[0]->xmax - array[0]->dx / 2.0 + dx / 2.0;
      ymax = array[0]->ymax - array[0]->dy / 2.0 + dy / 2.0;
      zmax = array[0]->zmax - array[0]->dz / 2.0 + dz / 2.0;
      array[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz);
    }
  return;
}

void MultiGrid::SetInitialConditions(Array3D* phi, Array3D* rho)
{
  int i, j, k, n, q, index, index2;
  int Channelkmax, ChannelStopkmax, FixedRegionkmax;
  Channelkmax = (int)(ChannelDepth / (PixelSize / (int)GridsPerPixel));
  ChannelStopkmax = (int)(ChannelStopDepth / (PixelSize / (int)GridsPerPixel));
  int PixX, PixY, FilledPixX, FilledPixY;
  double PixXmin, PixYmin;
  //double Vgradient;
  double ChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / pow(MICRON_PER_CM, 3);
  // ChargeFactor converts doping in cm^-3 into the appropriate units
  double GridSpacing = PixelSize / (double)GridsPerPixel;
  double ChannelStopCharge, ChannelCharge, CollectCharge, FixedRegionCharge;
  ChannelStopCharge = ChannelStopDoping * MICRON_PER_CM / ((double)ChannelStopkmax * GridSpacing) * ChargeFactor;
  ChannelCharge =  ChannelDoping * MICRON_PER_CM / ((double)Channelkmax * GridSpacing) * ChargeFactor;
  // Background charge:

  //printf("In SetICs, Channelkmax = %d, ChannelStopkmax = %d, BackgroundCharge = %.4f, ChannelStopCharge = %.4f, ChannelCharge = %.4f\n",Channelkmax, ChannelStopkmax, BackgroundDoping*ChargeFactor,ChannelStopCharge, ChannelCharge);
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
	{
	  for (k=1; k<rho->nz-1; k++)
	    {
	      index = i + j * rho->nx + k * rho->nx * rho->ny;
	      rho->data[index] = BackgroundDoping * ChargeFactor;
	    }
	}
    }
  printf("Finished setting background charge, \n");
  fflush(stdout);
  // Potential on top
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  index = i + j * phi->nx + (phi->nz - 1) * phi->nx * phi->ny;
	  phi->data[index] = Vbb;
	}
    }
    printf("Finished setting Vtop, \n");
    fflush(stdout);
    // Fixed Potentials on bottom
  for (n=0; n<NumberofFixedRegions; n++)
    {
      FixedRegionkmax = (int)(FixedRegionChargeDepth[n] / (PixelSize / (int)GridsPerPixel));
      FixedRegionCharge =  FixedRegionDoping[i] * MICRON_PER_CM / ((double)FixedRegionkmax * GridSpacing) * ChargeFactor;
      for (i=0; i<phi->nx; i++)
	{
	  for (j=0; j<phi->ny; j++)
	    {
	      index = i + j * phi->nx;
	      if (phi->x[i] >= FixedRegionLowerLeft[n][0] && phi->x[i] <= FixedRegionUpperRight[n][0] && phi->y[j] >= FixedRegionLowerLeft[n][1] && phi->y[j] <= FixedRegionUpperRight[n][1])
		{
		  phi->data[index] = FixedRegionVoltage[n];
		  for (k=1; k<FixedRegionkmax+1; k++)
		    {
		      index2 = index + k * phi->nx * phi->ny;
		      rho->data[index2] = FixedRegionCharge;
		    }

		}
	    }
	}
    }
    printf("Finished setting Fixed Potentials, \n");
    fflush(stdout);
    // Potentials and Charges in Pixel Region
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  if (phi->x[i] < PixelRegionLowerLeft[n][0] || phi->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<phi->ny; j++)
	    {
	      if (phi->y[j] < PixelRegionLowerLeft[n][1] || phi->y[j] > PixelRegionUpperRight[n][1])
		{ 
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * phi->nx;
	      PixX = (int)((phi->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
	      PixY = (int)((phi->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[n][0] + PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[n][1] + PixY * PixelSize;
	      // First set the charges
	      //if (i%17 ==0 && j%17 ==0){printf("i = %d, j = %d, x = %.2f, y = %.2f, PixX = %d, pixY = %d, PixXmin = %.2f, PixYmin = %.2f\n",i,j,phi->x[i],phi->y[j],PixX,PixY,PixXmin,PixYmin);}
	      if (phi->y[j] <= PixYmin + ChannelStopWidth/2.0 || phi->y[j] >= PixYmin + PixelSize - ChannelStopWidth/2.0)
		{
		  // This is the Channel Stop Region
		  for (k=1; k<ChannelStopkmax+1; k++)
		    {
		      index2 = index + k * phi->nx * phi->ny;
		      rho->data[index2] = ChannelStopCharge;
		    }
		}
	      else
		{
		  // This is the channel region
		  for (k=1; k<Channelkmax+1; k++)
		    {
		      index2 = index + k * phi->nx * phi->ny;
		      rho->data[index2] = ChannelCharge;
		    }
		  for (q=0; q<NumberofFilledWells[n]; q++)
		    {
		      if (DistributedCharge == 1)
			{
			  CollectCharge = -CollectedCharge[n][q] / (3.33 * (10.0 - ChannelStopWidth) * GridSpacing) * pow(MICRON_PER_CM, 3) * ChargeFactor;
			}
		      else
			{
			  CollectCharge = -CollectedCharge[n][q] / (2.0 * GridSpacing * 2.0 * GridSpacing * GridSpacing) * pow(MICRON_PER_CM, 3) * ChargeFactor;// For point charge
			}
			  FilledPixX = (int)((FilledPixelCoords[n][q][0] - PixelRegionLowerLeft[n][0]) / PixelSize);
		      FilledPixY = (int)((FilledPixelCoords[n][q][1] - PixelRegionLowerLeft[n][1]) / PixelSize);
		      if (!(FilledPixX == PixX && FilledPixY == PixY))
			{
			  continue;
			}
		      if (DistributedCharge == 1)
			{
			  if (phi->x[i] >= PixXmin + PixelSize/3.0 && phi->x[i] <= PixXmin + 2.0*PixelSize/3.0)
			    {
			      index2 = index + (Channelkmax + 1) * phi->nx * phi->ny;
			      rho->data[index2] = CollectCharge;
			    }
			}
		      else //For point charge
			{
			  if (phi->x[i] >= PixXmin + PixelSize/2.0 - GridSpacing - .0001 && phi->x[i] <= PixXmin + PixelSize/2.0 + GridSpacing +.0001 && phi->y[j] >= PixYmin + PixelSize/2.0 - GridSpacing - .0001 && phi->y[j] <= PixYmin + PixelSize/2.0 + GridSpacing +.0001) // For point charge
			    {
			      index2 = index + (Channelkmax + 1) * phi->nx * phi->ny;
			      rho->data[index2] = CollectCharge;
			    }
			}
		    }
		}
	      // Now the gate voltages
	      if (phi->x[i] >= PixXmin + PixelSize/3.0 && phi->x[i] <= PixXmin + 2.0*PixelSize/3.0)
		{
		  // This is the collection region
		  phi->data[index] = Vparallel_hi;
		}
	      else
		{
		  // This is the barrier gate region
		  phi->data[index] = Vparallel_lo;
		}
	      if ((phi->y[j] <= PixYmin + ChannelStopWidth/2.0 || phi->y[j] >= PixYmin + PixelSize - ChannelStopWidth/2.0) && UndepletedChannelStop == 1)
		{
		  // This is the Channel Stop Region
		  phi->data[index] = Vchannelstop;
		}
	    }
	}
    }

  printf("Finished setting Pixel Potentials, \n");
  fflush(stdout);
  /*
  // Set average potential as starting point (Boundary conditions)
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  for (k=1; k<phi->nz-1; k++)
	    {
	      index = i + j * phi->nx + k * phi->nx * phi->ny;
	      Vgradient = Vaverage + float(k)/float(phi->nz - 1) * (Vbb - Vaverage);
	      phi->data[index] = Vgradient;
	    }
	}
	}*/
  return;
}


void MultiGrid::SOR(Array3D* phi, Array3D* rho, double w)
{
  // Assumes fixed potentials on the top and bottom, free BC on the sides.
  double newphi;
  double omw, w6, hsquared;
  int i, j, k, im, ip, j0, jm, jp, nxy, ind, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  omw = 1.0 - w; w6 = w / 6.0;
  for (i=0; i<phi->nx; i++)
    {
      if (i == 0) {im = 0;} else {im = 1;}
      if (i == phi->nx-1) {ip = 0;} else {ip = 1;}
      for (j=0; j<phi->ny; j++)
	{
	  j0 = j * phi->nx;
	  if (j == 0) {jm = 0;} else {jm = phi->nx;}
	  if (j == phi->ny-1) {jp = 0;} else {jp = phi->nx;}

	  for (k=1; k<phi->nz-1; k++)
	    {
	      ind = i + j0 + k * nxy;
	      indmx = ind - im;
	      indpx = ind + ip;
	      indmy = ind - jm;
	      indpy = ind + jp;
	      indmz = ind - nxy;
	      indpz = ind + nxy;
	      newphi = omw * phi->data[ind] + w6 * (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->data[indmz]+phi->data[indpz] + hsquared * rho->data[ind]);
	      phi->data[ind] = newphi;
	    }
	}
    }
  return;
}

void MultiGrid::RedBlack(Array3D* phi, Array3D* rho)
{
  // Assumes fixed potentials on the top and bottom, free BC on the sides.
  double newphi;
  double hsquared;
  int i, j, k, pass, im, ip, j0, jm, jp, nxy, ind, indmx, indpx, indmy, indpy, indmz, indpz;
  int kstart;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  for (pass=0; pass<2; pass++)
    { 
      for (i=0; i<phi->nx; i++)
 	{ 
	  if (i == 0) {im = 0;} else {im = 1;}
  	  if (i == phi->nx-1) {ip = 0;} else {ip = 1;}
	  for (j=0; j<phi->ny; j++)
	    {
	      if ((i+j)%2 == pass)
		{ 
		  kstart = 2;
		} 
	      else 
		{ 
		  kstart = 1;
		  }
	      j0 = j * phi->nx;
	      if (j == 0) {jm = 0;} else {jm = phi->nx;}
	      if (j == phi->ny-1) {jp = 0;} else {jp = phi->nx;}
	      for (k=kstart; k<phi->nz-1; k+=2)
		{
		      ind = i + j0 + k * nxy;
		      indmx = ind - im;
		      indpx = ind + ip;
		      indmy = ind - jm;
		      indpy = ind + jp;
		      indmz = ind - nxy;
		      indpz = ind + nxy;
		      newphi = (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->data[indmz]+phi->data[indpz] + hsquared * rho->data[ind]) / 6.0;
		      phi->data[ind] = newphi;
		}
	    }
	}
    }
  //printf("Grid Size = %d, error = %.3g\n",phi->nx, error);
  return;
}

double MultiGrid::Error(Array3D* phi, Array3D* rho)
{
  // Assumes fixed potentials on the top and bottom, free BC on the sides.
  double newphi, error = 0.0;
  double hsquared;
  int i, j, k, nxy, ind, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  for (i=1; i<phi->nx-1; i++)
    { 
      for (j=1; j<phi->ny-1; j++)
	{
	  for (k=1; k<phi->nz-1; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      indmx = ind - 1;
	      indpx = ind + 1;
	      indmy = ind - phi->nx;
	      indpy = ind + phi->nx;
	      indmz = ind - nxy;
	      indpz = ind + nxy;
	      newphi = (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->data[indmz]+phi->data[indpz] + hsquared * rho->data[ind]) / 6.0;
	      error = max(error, fabs(phi->data[ind] - newphi));
	    }
	}
    }
  //printf("Grid Size = %d, error = %.3g\n",phi->nx, error);
  return error;
}

void MultiGrid::Restrict(Array3D* phi, Array3D* newphi, Array3D* rho, Array3D* newrho)
{
  // Assumes fixed potentials on the top and bottom, free BC on the sides.
  int i, j, k, im, i0, ip, jm, j0, jp, km, k0, kp, nxy, newnxy;
  int newindex, ind000;
  int ind00p, ind0p0, indp00;
  int ind00m, ind0m0, indm00;
  int ind0mm, ind0mp, ind0pm, ind0pp;
  int indm0m, indm0p, indp0m, indp0p;
  int indmm0, indmp0, indpm0, indpp0;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  double phisum, rhosum;
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;

  for (i=0; i<newphi->nx; i++)
    {
      if (i == 0) {i0 = 0; im = 0; ip = 1;}
      else if (i == newphi->nx - 1) {i0 = phi->nx - 1; im = phi->nx - 2; ip = phi->nx - 1;}
      else {i0 = 2 * i; im = 2 * i - 1; ip = 2 * i + 1;}
      for (j=0; j<newphi->ny; j++)
	{
	  if (j == 0) {j0 = 0; jm = 0; jp = 1;}
	  else if (j == newphi->ny - 1) {j0 = phi->ny - 1; jm = phi->ny - 2; jp = phi->ny - 1;}
	  else {j0 = 2 * j; jm = 2 * j - 1; jp = 2 * j + 1;}
	  for (k=0; k<newphi->nz; k++)
	    {
	      if (k == 0) {k0 = 0; km = 0; kp = 0;}
	      else if (k == newphi->nz - 1) {k0 = phi->nz - 1; km = phi->nz - 1; kp = phi->nz - 1;}
	      else {k0 = 2 * k; km = 2 * k - 1; kp = 2 * k + 1;}
	      newindex = i + j*newphi->nx + k*newnxy;
	      ind000 = i0 + j0 * phi->nx + k0 * nxy;
	      phisum = 8.0 * phi->data[ind000];
	      rhosum = 8.0 * rho->data[ind000];
	      indm00 = im + j0 * phi->nx + k0 * nxy;
	      indp00 = ip + j0 * phi->nx + k0 * nxy;
	      ind0m0 = i0 + jm * phi->nx + k0 * nxy;
	      ind0p0 = i0 + jp * phi->nx + k0 * nxy;
	      ind00m = i0 + j0 * phi->nx + km * nxy;
	      ind00p = i0 + j0 * phi->nx + kp * nxy;
	      phisum += 4.0 * (phi->data[indm00] + phi->data[indp00] + phi->data[ind0m0] + phi->data[ind0p0] + phi->data[ind00m] + phi->data[ind00p]);
	      rhosum += 4.0 * (rho->data[indm00] + rho->data[indp00] + rho->data[ind0m0] + rho->data[ind0p0] + rho->data[ind00m] + rho->data[ind00p]);
	      indmm0 = im + jm * phi->nx + k0 * nxy;
	      indmp0 = im + jp * phi->nx + k0 * nxy;
	      indpm0 = ip + jm * phi->nx + k0 * nxy;
	      indpp0 = ip + jp * phi->nx + k0 * nxy;
	      phisum += 2.0 * (phi->data[indmm0] + phi->data[indmp0] + phi->data[indpm0] + phi->data[indpp0]);
	      rhosum += 2.0 * (rho->data[indmm0] + rho->data[indmp0] + rho->data[indpm0] + rho->data[indpp0]);
	      indm0m = im + j0 * phi->nx + km * nxy;
	      indm0p = im + j0 * phi->nx + kp * nxy;
	      indp0m = ip + j0 * phi->nx + km * nxy;
	      indp0p = ip + j0 * phi->nx + kp * nxy;
	      phisum += 2.0 * (phi->data[indm0m] + phi->data[indm0p] + phi->data[indp0m] + phi->data[indp0p]);
	      rhosum += 2.0 * (rho->data[indm0m] + rho->data[indm0p] + rho->data[indp0m] + rho->data[indp0p]);
	      ind0mm = i0 + jm * phi->nx + km * nxy;
	      ind0mp = i0 + jm * phi->nx + kp * nxy;
	      ind0pm = i0 + jp * phi->nx + km * nxy;
	      ind0pp = i0 + jp * phi->nx + kp * nxy;
	      phisum += 2.0 * (phi->data[ind0mm] + phi->data[ind0mp] + phi->data[ind0pm] + phi->data[ind0pp]);
	      rhosum += 2.0 * (rho->data[ind0mm] + rho->data[ind0mp] + rho->data[ind0pm] + rho->data[ind0pp]);
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      phisum += phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp];
	      rhosum += rho->data[indmmm] + rho->data[indmmp] + rho->data[indmpm] + rho->data[indmpp] + rho->data[indpmm] + rho->data[indpmp] + rho->data[indppm] + rho->data[indppp];
	      newphi->data[newindex] = phisum / 64.0;
	      newrho->data[newindex] = rhosum / 64.0;
	    }
	}
    }
  return;
}


void MultiGrid::Prolongate(Array3D* phi, Array3D* newphi)
{
  // Assumes fixed potentials on the top and bottom, free BC on the sides.
  int i, j, k, im, ip, jm, jp, km, kp, nxy, newnxy;
  int newindex;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;
  for (i=0; i<newphi->nx; i++)
    {
      im = max(0,rint((float)i / 2.0 - 0.1));
      ip = min(newphi->nx-1,rint((float)i / 2.0 + 0.1));
      for (j=0; j<newphi->ny; j++)
	{
	  jm = max(0,rint((float)j / 2.0 - 0.1));
	  jp = min(newphi->ny-1,rint((float)j / 2.0 + 0.1));
	  for (k=1; k<newphi->nz-1; k++)
	    {
	      km = rint((float)k / 2.0 - 0.1);
	      kp = rint((float)k / 2.0 + 0.1);
	      newindex = i + j * newphi->nx + k * newnxy;
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      newphi->data[newindex] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	}
    }
  return;
}

void MultiGrid::VCycle(Array3D** phi, Array3D** rho, double w, int nsteps, int ncycle)
{
  // Iterates a few steps (4) at each grid on the way down to smooth, then
  // iterates the coarsest grid down to machine precision,
  // then does a given number of steps (ncycle) at each finer scale on the way up.
  int i, j, niter;
  double error;
  
  for (i=0; i<nsteps; i++)
    {
      Restrict(phi[i], phi[i+1], rho[i], rho[i+1]);
      for (j=0; j<4; j++)
	{
	  //RedBlack(phi[i], rho[i]);
	  SOR(phi[i+1], rho[i+1], w);
	}
    }
  error = 100.0; niter = 0;
  while (error > 1.0E-12)
    {
      niter += 1;
      //RedBlack(phi[nsteps], rho[nsteps]);
      SOR(phi[nsteps], rho[nsteps], w);
      error = Error(phi[nsteps], rho[nsteps]);
    }
  printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g\n",phi[nsteps]->nx-1,phi[nsteps]->ny-1,phi[nsteps]->nz-1,niter,error);
  fflush(stdout);
  for (i=nsteps; i>0; i--)
    {
      Prolongate(phi[i], phi[i-1]);
      niter = ncycle * (int)pow(2,i-1);
      for (j=0; j<niter; j++)
	{
	  //RedBlack(phi[i-1], rho[i-1]);
	  SOR(phi[i-1], rho[i-1], w);
	}
      error = Error(phi[i-1], rho[i-1]);
      printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g\n",phi[i-1]->nx-1,phi[i-1]->ny-1,phi[i-1]->nz-1,niter,error);
      fflush(stdout);
    }
  return;
}


void MultiGrid::WriteOutputFile(string outputfiledir, string filenamebase, string name, Array3D* array)
{
  // This writes the data to the HDF files
  string underscore = "_", slash = "/", hdfname, filename;
  int* int_attr_data  = new int[3];
  double* double_attr_data  = new double[3];
  double* flipped_data  = new double[array->nx*array->ny*array->nz];
  int i, j, k, index, flipped_index;
  
  // There must be a better way to change from x fast to z fast, but this works.
  for (i=0; i<array->nx; i++)
    {
      for (j=0; j<array->ny; j++)
	{
	  for (k=0; k<array->nz; k++)
	    {
	      index = i + j * array->nx + k * array->nx * array->ny;
	      flipped_index = k + j * array->nz + i * array->nz * array->ny;
	      flipped_data[flipped_index] = array->data[index];
	    }
	}
    }
  hdfname = filenamebase+underscore+name;
  filename = outputfiledir+slash+hdfname;
  WriteHDF5File3(filename, hdfname, array->nx, array->ny, array->nz, flipped_data);
  // Now we write the attributes
  int_attr_data[0] = array->nx; int_attr_data[1] = array->ny; int_attr_data[2] = array->nz;
  WriteHDF5IntAttribute(filename, hdfname, "Dimension", 3, int_attr_data);
  double_attr_data[0] = Xmin; double_attr_data[1] = Ymin; double_attr_data[2] = Zmin;
  WriteHDF5DoubleAttribute(filename, hdfname, "Lower_Left", 3, double_attr_data);
  double_attr_data[0] = Xmax; double_attr_data[1] = Ymax; double_attr_data[2] = Zmax;
  WriteHDF5DoubleAttribute(filename, hdfname, "Upper_Right", 3, double_attr_data);

  delete[] flipped_data;
  return;
}

void MultiGrid::Gradient(Array3D* phi, Array3D** E)
{
  int i, j, k, n, nxy, ind;
  nxy = phi->nx * phi->ny;
  for (n=0; n<3; n++)
    {
      E[n] = new Array3D(phi->xmin, phi->xmax, phi->nx, phi->ymin, phi->ymax, phi->ny, phi->zmin, phi->zmax, phi->nz);
    }

  // Ex 
  for (j=0; j<phi->ny; j++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind]) / phi->dx;
	  ind = 1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-2 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind] - phi->data[ind-1]) / phi->dx;
	  for (i=2; i<phi->nx-2; i++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[0]->data[ind] = (-phi->data[ind+2] + 8.0 * phi->data[ind+1] - 8.0 * phi->data[ind-1] + phi->data[ind-2]) / (12.0 * phi->dx);
	    }
	}
    }

    // Ey 
  for (i=0; i<phi->nx; i++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = i + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind]) / phi->dy;
	  ind = i + phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-2) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-1) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind] - phi->data[ind-phi->nx]) / phi->dy;
	  for (j=2; j<phi->ny-2; j++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[1]->data[ind] = (-phi->data[ind+2*phi->nx] + 8.0 * phi->data[ind+phi->nx] - 8.0 * phi->data[ind-phi->nx] + phi->data[ind-2*phi->nx]) / (12.0 * phi->dy);
	    }
	}
    }

  // Ez 
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  ind = i + j * phi->nx;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind]) / phi->dz;
	  ind = i + j * phi->nx + nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dz);
	  ind = i + j * phi->nx + (phi->nz-2) * nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dz);
	  ind = i + j * phi->nx + (phi->nz-1) * nxy;
	  E[2]->data[ind] = (phi->data[ind] - phi->data[ind-nxy]) / phi->dz;
	  for (k=2; k<phi->nz-2; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[2]->data[ind] = (-phi->data[ind+2*nxy] + 8.0 * phi->data[ind+nxy] - 8.0 * phi->data[ind-nxy] + phi->data[ind-2*nxy]) / (12.0 * phi->dz);
	    }
	}
    }
  return;
}

void MultiGrid::Trace(Array3D* phi, Array3D** E, double* point)
{
  int i, nsteps = 0;
  double norm_step, E2;
  double x, y, z;
  x = point[0]; y = point[1]; z = point[2];
  double*  E_interp = new double[3];
  while (point[2] > phi->z[2] and nsteps < 1000)
    {
      nsteps += 1;
      E2 = 0.0;
      for (i=0; i<3; i++)
	{
	  E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
	  E2 += E_interp[i] * E_interp[i];
	}
      norm_step = phi->dz / sqrt(E2);
      for (i=0; i<3; i++)
	{
	  point[i] += E_interp[i] * norm_step;
	}
    }
  delete[] E_interp;
  return;
}

void MultiGrid::TracePaths(Array3D* phi, Array3D** E, double* point, ofstream& file)
{
  int i, nsteps = 0;
  double norm_step, E2;
  double x, y, z;
  x = point[0]; y = point[1]; z = point[2];
  double*  E_interp = new double[3];
  while (point[2] > phi->z[2] and nsteps < 1000)
    {
      nsteps += 1;
      E2 = 0.0;
      for (i=0; i<3; i++)
	{
	  E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
	  E2 += E_interp[i] * E_interp[i];
	}
      norm_step = phi->dz / sqrt(E2);
      for (i=0; i<3; i++)
	{
	  point[i] += E_interp[i] * norm_step;
	}
      file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
    }
  delete[] E_interp;
  return;
}

void MultiGrid::TraceGrid(Array3D* phi, Array3D** E)
{
  double x, y, z;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string filename = (outputfiledir+slash+outputfilebase+underscore+name);
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  file  << setw(15) << "xin" << setw(15) << "yin" << setw(15) << "zin" << setw(15) << "xout" << setw(15) << "yout" << setw(15) << "zout" << endl;

  x = PixelBoundaryLowerLeft[0] + PixelBoundaryStepSize[0] / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      y = PixelBoundaryLowerLeft[1] + PixelBoundaryStepSize[1] / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  point[0] = x;
	  point[1] = y;
	  z = phi->z[phi->nz-1];
	  point[2] = z;
	  if (LogPixelPaths == 0)
	    {
	      Trace(phi, E, point);
	    }
	  else
	    {
	      TracePaths(phi, E, point, file);
	    }
	      file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
	  y += PixelBoundaryStepSize[1];
	}
      x += PixelBoundaryStepSize[0];
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  return;
}

void MultiGrid::AddDipolePotentials(Array3D* phi)
{
  int i, j, k, ndipole, nimage, index;
  double r2, r3, deltaphi, p, px, py, pz, dipolesign;
  double DipoleFactor =  (QE/(4.0 * pi * EPSILON_0*EPSILON_SI)) * MICRON_PER_M;
  // DipoleFactor converts (charge in e-) * (separation in microns) into the appropriate units
  for (i=0; i<phi->nx; i++)
    { 
      for (j=0; j<phi->ny; j++)
	{
	  for (k=1; k<phi->nz-1; k++)
	    {
	      index = i + j * phi->nx + k * phi->nx * phi->ny;
	      for (ndipole=0; ndipole<NumberofDipoles; ndipole++)
		{
		  p = DipoleFactor * (double)DipoleCharge[ndipole] * (2.0 * DipoleZLocation[ndipole]);
		  // p is the dipole strength
		  // p(x,y,z) are the vector coordinates from the dipole to the grid point 
		  px = phi->x[i] - DipoleCoords[ndipole][0];
		  py = phi->y[j] - DipoleCoords[ndipole][1];
		  r2 = px*px + py*py;
		  for (nimage= - NumberofDipoleImages; nimage<NumberofDipoleImages+1; nimage++)
		    {
		      if (nimage > 0){ dipolesign = 1.0; }
		      else { dipolesign = -1.0; }
		      pz = phi->z[k] - (phi->zmin + 2.0 * (phi->zmax - phi->zmin) * double(nimage));
		      r3 = pow(r2 + pz*pz, 1.5);
		      deltaphi = p * dipolesign * pz / r3;
		      /*
		      if (i == 128 && j == 128 && k == 5)
			{
			  printf("i = %d, j = %d, k = %d, ndipole = %d, nimage = %d, deltaphi = %g\n",i,j,k,ndipole,nimage,deltaphi);
			  }*/
		      phi->data[index] += deltaphi;
		    }
		}
	    }
	}
    }
  return;
}

