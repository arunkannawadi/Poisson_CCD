/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** poisson.cpp **************************

#include "poisson.h"

//***************MAIN PROGRAM***********************

int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      printf("\nwrong number of arguments\n");
      printf("Only argument should be name of .cfg file");
      exit(0);
    }

  MultiGrid* multi = new MultiGrid(argv[1]);

  delete multi;
  return 0;
}
//***************END MAIN PROGRAM***********************

