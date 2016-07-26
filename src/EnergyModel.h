#ifndef ENERGYMODEL_H
#define ENERGYMODEL_H

/*
Energy Model data. Not meant to be its own class, rather a way of
organizing the energy model related components of an RNA object. As
such, the functions defined here are all on RNA objects. 
 */
#include <stdlib.h>
#include <stdio.h>

// forward declaration of RNA, for the compiler
struct RNA;
typedef struct RNA RNA;

#define NUM_HEXALOOPS 4
#define NUM_TLOOPS 30
#define NUM_TRILOOPS 2

typedef struct {
  char loop[6];
  double energy;
} TetraLoop;

typedef struct {
  char loop[5];
  double energy;
} TriLoop;

typedef struct {
  char loop[8];
  double energy;
} HexaLoop; 

typedef struct EnergyModel {
  // TODO:
  // * coaxial?
  // * coaxstack?
  double multiA;
  double multiB;
  double multiC;
  double auPenalty[4][4];
  double dangle3[4][4][4];
  double dangle5[4][4][4];
  HexaLoop hexaloop[NUM_HEXALOOPS];
  // * int11
  // * int21
  // * int22
  double hairpinLoop[30]; // loop
  double internalLoop[30]; // loop
  double bulgeLoop[30]; // loop
  // * miscloop
  double stack[4][4][4][4]; // g_stack
  TetraLoop tloop[NUM_TLOOPS];
  TriLoop triloop[NUM_TRILOOPS];
  double tstack[4][4][4][4];
  double int11[4][4][4][4][4][4];
  // * tstackcoax
  // * tstackh
  // * tstacki
  // * tstacki1n
  // * tstacki23
  // * tstackm
  double* scale;
  double* bscale;
} EnergyModel;

int baseMap(char* c);
void loadStack(RNA* strand); 
void loadHairpin(RNA* strand);
void loadInternal(RNA* strand);
void loadBulge(RNA* strand);
void load3dangle(RNA* strand);
void load5dangle(RNA* strand);
void loadTStack(RNA* strand);
void loadTloop(RNA* strand);
void loadHexaloop(RNA* strand);
void loadTriloop(RNA* strand);
void loadInt11(RNA* strand);
void initializeEnergyModel(RNA* strand);
void freeEnergyModel(RNA* strand);

/* double TScale(double T, double dG, double dH) {  */
/*   if (dG==INFINITE_ENERGY) return INFINITE_ENERGY; //keep infinity==infinity */
/*   return dH - (short)floor((float)(dH-dG)*T/T37inK+0.5); */
/* } */

#endif
