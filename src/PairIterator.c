#include "PairIterator.h"
#include <stdlib.h>
#include <stdio.h>

PairIterator* allocatePairIterator(int ms) {
  PairIterator* ret = (PairIterator*) malloc(sizeof(PairIterator));
  ret->size = 0;
  ret->iterator = 0;
  ret->maxSize = ms;
  ret->pairs = (int*) malloc(ms * sizeof(int));
  return ret;
}

void freePairIterator(PairIterator* pairIterator) {
  free(pairIterator->pairs);
  free(pairIterator);
}

int hasNext(PairIterator* pairIterator) {
  return pairIterator->size - pairIterator->iterator;
}

int next(PairIterator* pairIterator) {
  if(pairIterator->iterator >= pairIterator->size) {
    printf("Out of bounds next(pairIterator), index: %d, size: %d\n", 
	   pairIterator->iterator,
	   pairIterator->size);
    exit(1);
  }
  pairIterator->iterator++;
  return pairIterator->pairs[pairIterator->iterator];
}

int start(PairIterator* pairIterator) {
  pairIterator->iterator = 0;
  return pairIterator->pairs[0];
}

void add(PairIterator* pairIterator, int element) {
  if(pairIterator->size >= pairIterator->maxSize) {
    printf("Out of bounds add(pairIterator, i), size: %d, maxSize: %d\n", 
	   pairIterator->size,
	   pairIterator->maxSize);
    exit(1);
  }
  pairIterator->pairs[pairIterator->size] = element;
  pairIterator->size++;
}
