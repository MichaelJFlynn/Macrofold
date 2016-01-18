#include "RNA.h"
#include "EnergyModel.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

  newStrand->temperature = 37; // celcius

  initializeEnergyModel(newStrand);
  //printf("%g\n", newStrand->energyModel->stack[0][1][2][3]);

  return newStrand;
}

void freeRNA(RNA* strand) {
  free(strand->sequence);
  free(strand->intSequence);
  freeEnergyModel(strand);
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
  
  if(fgets(seqBuffer, MAX_SEQUENCE_LENGTH, fp) || (length = strlen(seqBuffer)) > MAX_SEQUENCE_LENGTH) {
    printf("Input sequence is longer than MAX_SEQUENCE_LENGTH: 10000\n");
    exit(1);
  }

  fclose(fp);

  sequence = (char*) calloc(length + 1, sizeof(char));
  memcpy(sequence, seqBuffer, length);

  return allocateRNA(sequence);
}

