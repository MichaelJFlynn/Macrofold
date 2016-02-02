#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include <string.h>
#include <math.h>

#define SCALE_PARAMETER 1

// A = 0, C = 1, G = 2, U = 3
int baseMap(char* letter) {
  char c = *letter;

  if(c == 'A' || c == 'a') {
    return 0;
  } else if( c == 'C' || c == 'c') {
    return 1;
  } else if( c == 'G' || c == 'g') {
    return 2;
  } else if( c == 'U' || c == 'u') {
    return 3;
  } else {
    printf("baseMap error, invalid base letter: %c\n", c);
    exit(1);
  }
}

void initializeEnergyModel(RNA* strand) {
  int i, size;
  double* filler;

  // cheap, lazy, perhaps dangerous hack to initialize EnergyModel to
  // contain INFINITY's in all the missing spaces in the arrays.
  strand->energyModel = malloc(sizeof(EnergyModel));
  if(! (sizeof(EnergyModel) % sizeof(double)) == 0) {
    printf("EnergyModel is not an even multiple of doubles\n");
    exit(1);
  }
  size = sizeof(EnergyModel) / sizeof(double);
  filler = (double*) strand->energyModel;
  for(i = 0; i < size; i++) {
    filler[i] = INFINITY;
  }



  strand->energyModel->multiA = 1;
  strand->energyModel->multiB = 1;
  strand->energyModel->multiC = 1;
  
  strand->energyModel->scale = (double*) malloc((strand->length + 1) * sizeof(double));
  strand->energyModel->scale[0] = 1;
  for(i = 1; i <= strand->length; i++) {
    strand->energyModel->scale[i] = strand->energyModel->scale[i-1] * SCALE_PARAMETER;
  }

  strand->energyModel->bscale = (double*) malloc((strand->length + 1) * sizeof(double));
  strand->energyModel->bscale[0] = 1;
  for(i = 1; i <= strand->length; i++) {
    strand->energyModel->bscale[i] = strand->energyModel->bscale[i-1] * strand->energyModel->multiB / SCALE_PARAMETER;
  }

  loadStack(strand); 
  loadHairpin(strand);
  loadInternal(strand);
  loadBulge(strand);
  load3dangle(strand);
  load5dangle(strand);
  loadTStack(strand);
  loadTloop(strand);
  loadHexaloop(strand);
  loadTriloop(strand);
}

void freeEnergyModel(RNA* strand) {
  free(strand->energyModel->scale);
  free(strand->energyModel->bscale);
  free(strand->energyModel);
}


void loadStack(RNA* strand) { 
  DataFile* stackData = readCSV("../data/stack.csv");
  int i;

  for(i = 0; i < stackData->nrow; i++) {
      strand->energyModel->stack[baseMap(get(stackData,i,0))]
	[baseMap(get(stackData,i,1))]
	[baseMap(get(stackData,i,2))]
	[baseMap(get(stackData,i,3))] = atof(get(stackData,i,4));
  }  
  freeDataFile(stackData);
}

void loadHairpin(RNA* strand) {
  DataFile* hairpinData = readCSV("../data/hairpin.csv");
  int i;
  // hardcoded in length
  for(i = 0; i < hairpinData->nrow; i++) {
    strand->energyModel->hairpinLoop[atoi(get(hairpinData, i, 0))] = atof(get(hairpinData, i,1));
  }

  freeDataFile(hairpinData);
}

void loadInternal(RNA* strand) {
  DataFile* internalData = readCSV("../data/internal.csv");
  int i;
  // hardcoded in length
  for(i = 0; i < internalData->nrow; i++) {
    strand->energyModel->internalLoop[atoi(get(internalData, i, 0))] = atof(get(internalData, i, 1));
  }
  freeDataFile(internalData);
}

void loadBulge(RNA* strand) {
  DataFile* bulgeData = readCSV("../data/bulge.csv");
  int i;
  for(i = 0; i < bulgeData->nrow; i++) {
    strand->energyModel->bulgeLoop[atoi(get(bulgeData, i, 0))] = atof(get(bulgeData,i,1));
  }

  freeDataFile(bulgeData);
}

void load3dangle(RNA* strand) {
  DataFile* dangleData = readCSV("../data/dangle3.csv");
  int i;
  for(i=0; i < dangleData->nrow; i++) {
    strand->energyModel->dangle3
      [baseMap(get(dangleData, i, 0))]
      [baseMap(get(dangleData, i, 1))]
      [baseMap(get(dangleData, i, 2))] = atof(get(dangleData, i, 3));
  }
  freeDataFile(dangleData);
}

void load5dangle(RNA* strand) {
  DataFile* dangleData = readCSV("../data/dangle5.csv");
  int i;
  for(i=0; i < dangleData->nrow; i++) {
    strand->energyModel->dangle5
      [baseMap(get(dangleData, i, 0))]
      [baseMap(get(dangleData, i, 1))]
      [baseMap(get(dangleData, i, 2))] = atof(get(dangleData, i, 3));

  }
  freeDataFile(dangleData);
}

void loadTStack(RNA* strand) {
  DataFile* stackData = readCSV("../data/stack.csv");
  int i;

  for(i = 0; i < stackData->nrow; i++) {
    strand->energyModel->tstack
      [baseMap(get(stackData, i, 0))]
      [baseMap(get(stackData, i, 1))]
      [baseMap(get(stackData, i, 2))]
      [baseMap(get(stackData, i, 3))] = atof(get(stackData, i, 4));
  }
  freeDataFile(stackData);
}

void loadTloop(RNA* strand) {
  DataFile* tloopData = readCSV("../data/tloop.csv");
  int i;
  
  if(tloopData->nrow != NUM_TLOOPS) {
    printf("Number of tloop data entries does not match NUM_TLOOPS.\n");
    exit(1);
  }

  for(i = 0; i < tloopData->nrow; i++) {
    memcpy(strand->energyModel->tloop[i].loop, get(tloopData,i, 0), 6);
    strand->energyModel->tloop[i].energy = atof(get(tloopData, i, 1));
  }
  freeDataFile(tloopData);
}

void loadHexaloop(RNA* strand) {
  DataFile* hexaData = readCSV("../data/hexaloop.csv");
  int i;

  if(hexaData->nrow != NUM_HEXALOOPS) {
    printf("Number of hexaloop data entries does not match NUM_HEXALOOPS.\n");
    exit(1);
  }

  for(i = 0; i < hexaData->nrow; i++) {
    memcpy(strand->energyModel->hexaloop[i].loop, get(hexaData, i, 0), 8);
    strand->energyModel->hexaloop[i].energy =  atof(get(hexaData, i,1));
  }
  freeDataFile(hexaData);
}

void loadTriloop(RNA* strand) {
  DataFile* triData = readCSV("../data/triloop.csv");
  int i;

  if(triData->nrow != NUM_TRILOOPS) {
    printf("Number of triloop data entries does not match NUM_TRILOOPS.\n");
    exit(1);
  }

  for(i = 0; i < triData->nrow; i++) {
    memcpy(strand->energyModel->triloop[i].loop, get(triData, i, 0), 5);
    strand->energyModel->triloop[i].energy = atof(get(triData, i, 1));
  }
}

