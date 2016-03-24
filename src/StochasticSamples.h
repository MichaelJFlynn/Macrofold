
struct RNA;
typedef struct RNA RNA;

typedef struct StochasticSamples {

} StochasticSamples;

typedef struct stackNode {
  int i;
  int j;
  int matrix; /* [0, 1, 2] ~ [Q', Q1, Q] */
  struct stackNode* next;
} stackNode;

void sample(RNA* strand, int numSamples);
void push(stackNode** s, int i, int j, int m) {

}
