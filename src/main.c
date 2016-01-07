// Written by Mike Flynn

#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
  RNA* rna = allocateRNA();
  
  printf("RNA loaded.\n");
  
  // TODO: might have to change this to something different for every
  // operating system... damn
  DataFile* data = readCSV("..\\data\\stack.csv");
  freeDataFile(data);

  freeRNA(rna);
  printf("RNA freed.\n");
  return 0;
}
