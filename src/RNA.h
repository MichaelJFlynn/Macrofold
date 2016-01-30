#ifndef RNA_H
#define RNA_H

// forward declaration, for compiler
struct EnergyModel;
typedef struct EnergyModel EnergyModel;

struct AllowedPairs;
typedef struct AllowedPairs AllowedPairs;

typedef double** PartitionFunctionMatrix;
typedef double* PartitionFunctionVector;

typedef struct {
  double* Z;
  double** Zb;
  double** Z1;
  double** Z2;
} PartitionFunctionData;

typedef struct RNA {
  int length;
  char* sequence;
  // sequence mapped to integers (see EnergyModel.h -> baseMap):
  int* intSequence;
  double temperature;
  EnergyModel* energyModel;
  PartitionFunctionData pfData;
  AllowedPairs* allowedPairs;
} RNA;

RNA* allocateRNA(char* sequence);
RNA* readSequenceFile(char* filename);
void freeRNA(RNA* strand);
#endif
