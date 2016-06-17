#include "RNA.h"
#include "AllowedPairs.h"
#include "PairIterator.h"
#include "EnergyFunctions.h" // for TURN definition
#include "PartitionFunction.h"
#include <stdlib.h>
#include <stdio.h>

AllowedPairs* fromAllPairs(RNA* strand) {
  int length = strand->length;
  int i,j;
  AllowedPairs* ap = (AllowedPairs*) malloc(sizeof(AllowedPairs));
  ap->size = length;
  ap->thresh = 0;
  ap->ji = (PairIterator**) malloc(length * sizeof(PairIterator*));
  ap->ij = (PairIterator**) malloc(length * sizeof(PairIterator*));

  for(i =0; i < length; i++) {
    ap->ij[i] = allocatePairIterator(length);
    ap->ji[i] = allocatePairIterator(length);
  }

  for(i = 0; i < length;  i++) {
    for(j =i+TURN; j < length; j++) {
      if(isCannonical(strand, i, j)) {
	add(ap->ij[i], j);
      }
    }
  }

  for(j = length - 1; j >= 0;  j--) {
    for(i =j-TURN; i >= 0; i--) {
      if(isCannonical(strand, i, j)) {
	add(ap->ji[j], i);
      }
    }
  }
  return ap;
}


AllowedPairs* fromProbablePairs(RNA* strand, float threshold) { 
  int length = strand->length;
  int i,j;
  AllowedPairs* ap = (AllowedPairs*) malloc(sizeof(AllowedPairs));
  ap->size = length;
  ap->thresh = threshold;
  ap->ji = (PairIterator**) malloc(length * sizeof(PairIterator*));
  ap->ij = (PairIterator**) malloc(length * sizeof(PairIterator*));

  double** P = strand->partitionFunction->P;
  double** Zb = strand->partitionFunction->Zb;

  for(i =0; i < length; i++) {
    ap->ij[i] = allocatePairIterator(length);
    ap->ji[i] = allocatePairIterator(length);
  }

  for(i = 0; i < length;  i++) {
    for(j =i+TURN; j < length; j++) {
      if(isCannonical(strand, i, j) && P[i][j] > threshold) {
	add(ap->ij[i], j);
      }
    }
  }

  for(j = length - 1; j >= 0;  j--) {
    for(i =j-TURN; i >= 0; i--) {
      if(isCannonical(strand, i, j) && P[i][j] > threshold) {
	add(ap->ji[j], i);
      }
    }
  }
  return ap;  
}


AllowedPairs* fromProbablePairs2(RNA* strand, float threshold) { 
  int length = strand->length;
  int i,j;
  AllowedPairs* ap = (AllowedPairs*) malloc(sizeof(AllowedPairs));
  ap->size = length;
  ap->thresh = threshold;
  ap->ji = (PairIterator**) malloc(length * sizeof(PairIterator*));
  ap->ij = (PairIterator**) malloc(length * sizeof(PairIterator*));

  double** P = strand->partitionFunction->P;
  double** Zb = strand->partitionFunction->Zb;

  for(i =0; i < length; i++) {
    ap->ij[i] = allocatePairIterator(length);
    ap->ji[i] = allocatePairIterator(length);
  }

  for(i = 0; i < length;  i++) {
    for(j =i+TURN; j < length; j++) {
      if(isCannonical(strand, i, j) && Zb[i][j] > threshold) {
	add(ap->ij[i], j);
      }
    }
  }

  for(j = length - 1; j >= 0;  j--) {
    for(i =j-TURN; i >= 0; i--) {
      if(isCannonical(strand, i, j) && Zb[i][j] > threshold) {
	add(ap->ji[j], i);
      }
    }
  }
  return ap;  
}


void freeAllowedPairs(AllowedPairs* allowedPairs) {
  int i;
  for(i = 0; i < allowedPairs->size; i++) { 
    freePairIterator(allowedPairs->ij[i]);
    freePairIterator(allowedPairs->ji[i]);
  }
  free(allowedPairs);
}

void printAllowedPairs(AllowedPairs* ap) {
  int i, j, k;
  PairIterator* p;

  printf("ij ##################################\n");
  for(i = 0; i < ap->size; i++) {
    p = ap->ij[i];
    printf("%d: ", i);
    for(k = start(p); hasNext(p); k = next(p)) {
      printf("%d ", k);
    }
    printf("\n"); 
  }

  printf("ji ##################################\n");
  for(j = ap->size - 1; j >= 0; j--) {
    p = ap->ji[j];
    printf("%d: ", j);
    for(k = start(p); hasNext(p); k = next(p)) {
      printf("%d ", k);
    }
    printf("\n"); 
  }
}
