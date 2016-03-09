#ifndef PARTITION_FUNCTION_H
#define PARTITION_FUNCTION_H
// forward declaration of RNA, for the compiler
struct RNA;
typedef struct RNA RNA;

typedef struct PartitionFunction {
  double** Z;
  double** Zb;
  double** Z1;
  double** Z2;
  double** P;
  int filledZbZ1Z2;
  int filledZbZ1Z2extended;
  int filledZ;
  int filledP;
} PartitionFunction;

PartitionFunction* allocatePartitionFunction(int length);
void freePartitionFunction(RNA* strand);
void fillZbZ1Z2(RNA* strand);
void fillZ(RNA* strand); 
void fillExtendedZbZ1Z2(RNA* strand);
void fillP(RNA* strand);
#endif
