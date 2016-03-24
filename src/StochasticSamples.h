
struct RNA;
typedef struct RNA RNA;

typedef int* Structure;

typedef struct StochasticSamples {
  Structure* structures;
  int numSamples;
} StochasticSamples;

typedef struct stackNode {
  int i;
  int j;
  int matrix; /* [0, 1, 2] ~ [Q', Q1, Q] */
  struct stackNode* next;
} stackNode;


StochasticSamples* allocateStochasticSamples(RNA* strand, int numSamples);
void freeSamples(StochasticSamples* samples);

Structure allocateStructure(RNA* strand);
void freeStructure(Structure str);
void pair(Structure s, int i, int j);


void sample(RNA* strand, int numSamples);
void push(stackNode** s, int i, int j, int m);
