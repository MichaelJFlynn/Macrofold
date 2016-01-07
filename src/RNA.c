#include "RNA.h"
#include "EnergyModel.h"
#include <stdlib.h>

RNA* allocateRNA() {
  RNA* newStrand = (RNA*) malloc(sizeof(RNA));
  
  newStrand->length = 0;
  newStrand->temperature = 37; // celcius

  initializeEnergyModel(newStrand);
  //printf("%g\n", newStrand->energyModel->stack[0][1][2][3]);

  return newStrand;
}

void freeRNA(RNA* strand) {
  free(strand->energyModel);
  free(strand);
}

int initialized(RNA* strand) {
  return strand->length;
}

