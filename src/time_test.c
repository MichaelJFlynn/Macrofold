// Written by Mike Flynn 2016-06-02;
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
#include <time.h>

/* 
This program runs tests to make sure the build of Macrofold is
correct.
 */
int main(int argc, char* argv[]) {
  if(argc != 2) {
    printf("Usage:\ntime_test <seq_file>\n");
    return 0;
  }

  clock_t start, normalTime, fourTime,
    sixTime, eightTime, tenTime;

  RNA* strand = readSequenceFile(argv[1]);
  computePartitionFunction(strand);

  start = clock();
  sample(strand, 1000);
  normalTime = clock() - start;
  
  strand->allowedPairs = fromProbablePairs(strand, 1e-10);
  start = clock();
  sample(strand, 1000);
  tenTime = clock() - start;

  strand->allowedPairs = fromProbablePairs(strand, 1e-8);
  start = clock();
  sample(strand, 1000);
  eightTime = clock() - start;

  strand->allowedPairs = fromProbablePairs(strand, 1e-6);
  start = clock();
  sample(strand, 1000);
  sixTime = clock() - start;

  /* strand->allowedPairs = fromProbablePairs(strand, 1e-4); */
  /* start = clock(); */
  /* sample(strand, 1000); */
  /* fourTime = clock() - start; */

  printf("%s\t%d\t%ju\t%ju\t%ju\t%ju\n", 
	 argv[1], 
	 strand->length, 
	 normalTime, 
	 tenTime,
	 eightTime,
	 sixTime
	 );
}
