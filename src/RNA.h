#ifndef RNA_H
#define RNA_H

// forward declarations, for compiler
struct EnergyModel;
typedef struct EnergyModel EnergyModel;

struct AllowedPairs;
typedef struct AllowedPairs AllowedPairs;

struct PartitionFunction;
typedef struct PartitionFunction PartitionFunction;

typedef struct RNA {
  int length;
  char* sequence;
  // sequence mapped to integers (see EnergyModel.h -> baseMap):
  int* intSequence;
  double temperature;
  EnergyModel* energyModel;
  PartitionFunction* partitionFunction;
  AllowedPairs* allowedPairs;
} RNA;

RNA* allocateRNA(char* sequence);
RNA* readSequenceFile(char* filename);
void freeRNA(RNA* strand);
double getFreeEnergy(RNA* strand);
int isCannonical(RNA* strand, int i, int j);
#endif
