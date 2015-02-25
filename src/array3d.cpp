/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** array3d.cpp **************************

#include "array3d.h"

Array3D::Array3D(double Xmin, double Xmax, int Nx, double Ymin, double Ymax, int Ny, double Zmin, double Zmax, int Nz)
{
  int i, j, k;
  nx = Nx; ny = Ny; nz = Nz; 
  xmin = Xmin; xmax = Xmax; ymin = Ymin; ymax = Ymax; zmin = Zmin; zmax = Zmax;
  dx = (xmax - xmin) / (double) nx;
  dy = (ymax - ymin) / (double) ny;
  dz = (zmax - zmin) / (double) nz;
  volume = dx * dy * dz;
  x = new double[nx]; y = new double[ny]; z = new double[nz];
  data = new double[nx * ny * nz];
  for (k=0; k<nz; k++)
    {
      z[k] = zmin + dz/2.0 + (double) k * dz;
    }
  for (j=0; j<ny; j++)
    {
      y[j] = ymin + dy/2.0 + (double) j * dy;
    }
  for (i=0; i<nx; i++)
    {
      x[i] = xmin + dx/2.0 + (double) i * dx;	  
      for (j=0; j<ny; j++)
	{
	  for (k=0; k<nz; k++)
	    {
	      data[i + j*nx + k*nx*ny] = 0.0;
	    }
	}
    }
}

Array3D::~Array3D()
{
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] data;
}

double Array3D::PyramidalKernel3D(double deltax, double deltay, double deltaz)
{
  if (deltax>=1.0 || deltay>=1.0 || deltaz>=1.0)
    {
      return 0.0;
    }
  else
    {
      return (1.0-deltax)*(1.0-deltay)*(1.0-deltaz);
    }
}

double Array3D::DataInterpolate3D(double xprime, double yprime, double zprime)
{
  int i, j, k, m, n, ml, nl, p, pl;
  double d, norm, PK, deltax, deltay, deltaz; 
  i=(int)ceil((xprime-xmin)/dx);
  j=(int)ceil((yprime-ymin)/dy);
  k=(int)ceil((zprime-zmin)/dz);
  d=0.0; norm = 0.0;
  for (m=i-1; m<i+2; m++)
    {
      ml = max(0,min(nx-1,m));
      deltax=fabs((xprime-x[ml])/dx);
      for (n=j-1; n<j+2; n++)
	{
	  nl = max(0,min(ny-1,n));
	  deltay=fabs((yprime-y[nl])/dy);
	  for (p=k-1; p<k+2; p++)
	    {
	      pl = max(0,min(nz-1,p));
	      deltaz=fabs((zprime-z[pl])/dz);
	      PK = PyramidalKernel3D(deltax,deltay,deltaz);
	      norm += PK;
	      d += PK * data[ml + nl * nx + pl * nx * ny];
	    }
	}
    }
  if (norm < 1.0E-15)
    {
      return 0.0;
    }
  else
    {
      return d/norm;
    }
}

