/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/
//****************** fileio.h **************************

#include <stdio.h> 
#include <string.h>
#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace std;
double GetDoubleParam(string, string, double);
double* GetDoubleList(string, string, int, double*);
int GetIntParam(string, string, int);
int* GetIntList(string, string, int, int*);
string GetStringParam(string, string, string);
string ReadPar(string, string); 
