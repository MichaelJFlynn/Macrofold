// Written by Mike Flynn 2016-06-02
#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include "EnergyFunctions.h"
#include "PartitionFunction.h"
#include "AllowedPairs.h"
#include "MacrofoldConsole.h"
#include "StochasticSamples.h"
#include <stdio.h>
#include <math.h>

/* 
This program runs tests to make sure the build of Macrofold is
correct.
 */
int main(int argc, char* argv[]) {
  if(argc != 2) {
    printf("Usage:\ntest_pij <seq_file>\n");
    return 0;
  }

  RNA* strand = readSequenceFile(argv[1]);
  computePartitionFunction(strand);
  
  int i, j;
  for(i = 0; i < strand->length; i++) {
    //    for(j = 0; j < 2*strand->length; j++) {
    //if(i < j) {
    printf("%.4e\t", strand->partitionFunction->Z[0][i]);
    //} else {
    ///printf("           \t");
    //}
    //}
    //printf("\n");
  } printf("\n");
  
}
