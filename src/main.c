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
  fillExtendedZbZ1Z2(strand);
  fillP(strand);


  printf("%g\n", getFreeEnergy(strand));
  int i,j;
  for(i = 0; i < strand->length; i++) {
    for(j = 0; j < strand->length; j++) {
      double prob = strand->partitionFunction->P[i][j];
      if(prob > .001) {
	printf("P[%d][%d] = %g\n", i, j, prob);
      }
    }
  } 

  /*
    Reporting Free Energies
    Finding Pijs
    Stochastic traceback
   */

  freeRNA(strand);
  printf("RNA freed.\n");
  return 0;
}
