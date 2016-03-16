#ifndef ALLOWEDPAIRS_H
#define ALLOWEDPAIRS_H

struct PairIterator;
typedef struct PairIterator PairIterator;

struct RNA;
typedef struct RNA RNA;

typedef struct AllowedPairs {
  PairIterator** ij;
  PairIterator** ji;
  int size;
  float thresh;
} AllowedPairs;

void freeAllowedPairs(AllowedPairs* allowPairs);
AllowedPairs* fromAllPairs(RNA* strand);
void printAllowedPairs(AllowedPairs* ap);

#endif
