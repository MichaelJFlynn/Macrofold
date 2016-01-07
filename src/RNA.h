#ifndef RNA_H
#define RNA_H

// forward declaration, for compiler
struct EnergyModel;
typedef struct EnergyModel EnergyModel;

typedef double** PartitionFunctionMatrix;
typedef double* PartitionFunctionVector;

typedef struct {
  PartitionFunctionVector Z;
  PartitionFunctionMatrix Zb;
  PartitionFunctionMatrix Z1;
  PartitionFunctionMatrix Z2;
} PartitionFunctionData;

typedef struct RNA {
  int length;
  double temperature;
  EnergyModel* energyModel;
  PartitionFunctionData pfData;
} RNA;

int initialized(RNA* strand);

RNA* allocateRNA();
void freeRNA(RNA* strand);
#endif
