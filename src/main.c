// Written by Mike Flynn

#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include "EnergyFunctions.h"
#include <stdio.h>




int main(int argc, char* argv[]) {
  // tests 
  RNA* strand = readSequenceFile("sample.seq");
  
  printf("%g\n", hairpinTerm(strand, 3, 7));
  

  freeRNA(strand);
  printf("RNA freed.\n");
  return 0;
}
