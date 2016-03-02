// Written by Mike Flynn

#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include "EnergyFunctions.h"
#include "PartitionFunction.h"
#include <stdio.h>
#include <math.h>


int main(int argc, char* argv[]) {
  // tests 
  RNA* strand = readSequenceFile("sample.seq");
  fillZbZ1Z2(strand);
  fillZ(strand);
  

  printf("%g\n", getFreeEnergy(strand));


  /*
    Reporting Free Energies
    Finding Pijs
    Stochastic traceback
   */

  freeRNA(strand);
  printf("RNA freed.\n");
  return 0;
}
