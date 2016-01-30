
struct PairIterator;
typedef struct PairIterator PairIterator;

typedef struct AllowedPairs {
  PairIterator** ij;
  PairIterator** ji;
  int size;
  float thresh;
} AllowedPairs;

void freeAllowedPairs(AllowedPairs* allowPairs);
AllowedPairs* fromAllPairs(int length);
