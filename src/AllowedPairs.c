#include "AllowedPairs.h"
#include "PairIterator.h"
#include "EnergyFunctions.h" // for TURN definition
#include <stdlib.h>

AllowedPairs* fromAllPairs(int length) {
  int i,j;
  AllowedPairs* ap = (AllowedPairs*) malloc(sizeof(AllowedPairs));
  ap->size = length;
  ap->thresh = 0;
  ap->ji = (PairIterator**) malloc(length * sizeof(PairIterator*));
  ap->ij = (PairIterator**) malloc(length * sizeof(PairIterator*));

  for(i =0; i <= length; i++) {
    ap->ij[i] = allocatePairIterator(length);
    ap->ji[i] = allocatePairIterator(length);
  }

  for(i = 0; i < length;  i++) {
    for(j =i+TURN; j < length; j++) {
      add(ap->ij[i], j);
    }
  }

  for(j = length - 1; j >= 0;  j--) {
    for(i =j-TURN; i >= 0; i--) {
      add(ap->ji[j], i);
    }
  }
  return ap;
}

void freeAllowedPairs(AllowedPairs* allowedPairs) {
  freePairIterator(allowedPairs->ij);
  freePairIterator(allowedPairs->ji);
  free(allowedPairs);
}
