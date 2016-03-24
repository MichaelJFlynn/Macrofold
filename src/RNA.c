#include "RNA.h"
#include "EnergyModel.h"
#include "AllowedPairs.h"
#include "PartitionFunction.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_SEQUENCE_LENGTH 10000

RNA* allocateRNA(char* sequence) {
  RNA* newStrand = (RNA*) malloc(sizeof(RNA));
  int i;

  newStrand->sequence = sequence;
  newStrand->length = strlen(sequence);
  newStrand->intSequence = malloc(newStrand->length * sizeof(int));
  for(i=0; i < newStrand->length; i++) {
    newStrand->intSequence[i] = baseMap(&newStrand->sequence[i]);
  }

  newStrand->temperature = 37 + 273; // kelvin
  newStrand->allowedPairs = fromAllPairs(newStrand);
  initializeEnergyModel(newStrand);
  newStrand->partitionFunction = allocatePartitionFunction(newStrand->length);
  newStrand->samples = 0;
  //printf("%g\n", newStrand->energyModel->stack[0][1][2][3]);
  return newStrand;
}

void freeRNA(RNA* strand) {
  free(strand->sequence);
  free(strand->intSequence);
  freeEnergyModel(strand);
  freeAllowedPairs(strand->allowedPairs);
  freePartitionFunction(strand);
  free(strand);
}

RNA* readSequenceFile(char* filename) { 
  FILE* fp = fopen(filename, "r");
  char* sequence;
  int length;

  if(!fp) {
    printf("Sequence file %s not found.\n", filename);
    exit(1);
  }

  char seqBuffer[MAX_SEQUENCE_LENGTH];
  fgets(seqBuffer, MAX_SEQUENCE_LENGTH, fp);
  if( (length = strlen(seqBuffer)) > MAX_SEQUENCE_LENGTH) {
    printf("Input sequence is longer than MAX_SEQUENCE_LENGTH: %d\n", MAX_SEQUENCE_LENGTH);
    exit(1);
  }

  fclose(fp);

  sequence = (char*) calloc(length + 1, sizeof(char));
  memcpy(sequence, seqBuffer, length);
  return allocateRNA(sequence);
}

int isCannonical(RNA* strand, int i, int j) { 
  if(i < 0 || j < 0 || i >= strand->length || j >= strand->length) {
    printf("isCannonical error: out of bounds (%d,%d)\n", i, j);
    exit(1);
  }
  char ic = strand->sequence[i];
  char jc = strand->sequence[j];
  if(ic == 'A' && jc == 'U') {
    return 1;
  } else if(ic == 'U' && jc == 'A') {
    return 1;
  } else if(ic == 'G' && jc == 'C') {
    return 1;
  } else if(ic == 'C' && jc == 'G') {
    return 1;
  } else if(ic == 'G' && jc == 'U') {
    return 1;
  } else if(ic == 'U' && jc == 'G') {
    return 1;
  } else {
    return 0;
  }
}

/* Outputs the change in energy from the completely unfolded state to the folded state

 */
double getFreeEnergy(RNA* strand) {
  /* fputs("#T\t-RT ln Z\tZ\n", dGFile); */
  /* fprintf(dGFile, "%g\t%g\t%g\n", t, -RT * log(Q5(g_len) - Z0) - RT * g_len * log(g_scale), (Q5(g_len) - Z0) * g_scalen[g_len]); */
  double RT = .0019872 * strand->temperature;
  return -RT * log(strand->partitionFunction->Z[0][strand->length - 1]) - RT * strand->length * log(strand->energyModel->scale[1]);
}


