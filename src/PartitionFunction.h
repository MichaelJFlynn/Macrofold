#ifndef PARTITION_FUNCTION_H
#define PARTITION_FUNCTION_H
// forward declaration of RNA, for the compiler
struct RNA;
typedef struct RNA RNA;

typedef struct PartitionFunction {
  double* Z;
  double** Zb;
  double** Z1;
  double** Z2;
} PartitionFunction;

PartitionFunction* allocatePartitionFunction(int length);
void freePartitionFunction(RNA* strand);
void fillZbZ1Z2(RNA* strand);
void fillZ(RNA* strand); 

#endif
