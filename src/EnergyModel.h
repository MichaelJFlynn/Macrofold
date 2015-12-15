#ifndef ENERGYMODEL_H
#define ENERGYMODEL_H

#include <stdlib.h>

typedef struct {
  // TODO:
  // * coaxial?
  // * coaxstack?
  // * dangle
  // * hexaloop
  // * int11
  // * int21
  // * int22
  double hairpinLoop[30]; // loop
  double internalLoop[30]; // loop
  double bulgeLoop[30]; // loop
  // * miscloop
  double stack[4][4][4][4]; // g_stack
  // * tloop
  // * triloop
  // * tstack
  // * tstackcoax
  // * tstackh
  // * tstacki
  // * tstacki1n
  // * tstacki23
  // * tstackm
  double scale; 
  double scalePower[30];
  

} EnergyModel;


// A = 0, C = 1, G = 2, U = 3 
int baseMap(char c) {
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

RNA loadDataTables(RNA rna) {
  EnergyModel model = rna.energyModel;
  model->scale = 0.34; // kcal/mol/base 

  model->scalePower = (double*) malloc(sizeof(double), length);
  scalePower[0] = model->scale;

  char** bulgeData = read.csv("data/bulge.csv");
  char** hairpinData = read.csv("data/hairpin.csv");
  char** hexaloopData = read.csv("data/hexaloop.csv");
  char** internalData = read.csv("data/internal.csv");
  char** tloopData = read.csv("data/tloop.csv");
  char** triloopData = read.csv("data/triloop.csv");
  char** tstackData = read.csv("data/tstack.csv");
    
   


  for(int i = 1; i < length; i++) {
    scalePower

  }
  
  
}

void loadStack(RNA strand) {
  char** stackData = read.csv("data/stack.csv");
  
  if(!strand.energyModel) {
    printf("Energy model not initialized in loadStack");
    exit(1);
  }

  for(int i = 0; i < nrow(stackData); i++) {
    strand.energyModel->stack[baseMap(stackData[i][0]),
			      baseMap(stackData[i][1]),
			      baseMap(stackData[i][2]),
			      baseMap(stackData[i][3])
			      ] = atof(stackData[i][4]);
  }
}

void loadHairpin(RNA strand) {
  char** hairpinData = read.csv("data/hairpin.csv");
  
  // hardcoded in length
  for(int i = 0; i < 30; i++) {
    strand.energyModel->hairpinLoop[i] = atod(hairpinData[i][1]); 
  }
}

void loadInternal(RNA strand) {
  char** internalData = read.csv("data/internal.csv");
  
  // hardcoded in length
  for(int i = 0; i < 30; i++) {
    strand.energyModel->internalLoop[i] = atod(internalData[i][1]); 
  }
}

void loadBulge(RNA strand) {
  char** bulgeData = read.csv("data/bulge.csv");

  for(int i = 0; i < 30; i++) {
    strand.energyModel->bulgeLoop[i] = atod(internalData[i][1]);
  }
}


double TScale(double T, double dG, double dH) { 
  if (dG==INFINITE_ENERGY) return INFINITE_ENERGY; //keep infinity==infinity
  return dH - (short)floor((float)(dH-dG)*T/T37inK+0.5);
}

#endif
