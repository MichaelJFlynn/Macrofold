typedef struct PairIterator {
  int* pairs;
  int size;
  int iterator;
  int maxSize;
} PairIterator;

PairIterator* allocatePairIterator(int ms);
void freePairIterator();
int hasNext(PairIterator* pairIterator);
int next(PairIterator* pairIterator);
void start(PairIterator* pairIterator);
void add(PairIterator* pairIterator, int element);
