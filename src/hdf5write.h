/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** hdf5write.h **************************

#include <string>

#include "H5Cpp.h"
using namespace std;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

int WriteHDF5File3(string, string, int, int, int, double*); // Rank 3 files 
int WriteHDF5File2(string, string, int, int, double*); // Rank 2 files 
int WriteHDF5File1(string, string, int, double*); // Rank 1 files
int WriteHDF5IntAttribute(string, string, string, int, int*);
int WriteHDF5DoubleAttribute(string, string, string, int, double*);

