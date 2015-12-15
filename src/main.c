// Written by Mike Flynn

#include "loadData.h"
#include "RNA.h"

int main(int argc, char* argv[]) {
  
  RNA rna = allocateRNA();
  loadData();
  computePartitionFunction(rna);
  rna.Z1(1, 

}
