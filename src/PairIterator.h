#ifndef PAIRITERATOR_H
#define PAIRITERATOR_H

typedef struct PairIterator {
  int* pairs;
  int size;
  int iterator; // TODO: int* ?
  int maxSize;
} PairIterator;

PairIterator* allocatePairIterator(int ms);
void freePairIterator();
int hasNext(PairIterator* pairIterator);
int next(PairIterator* pairIterator);
int start(PairIterator* pairIterator);
void add(PairIterator* pairIterator, int element);
void printPairIterator(PairIterator* pairIterator);
#endif
