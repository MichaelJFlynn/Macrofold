#ifndef RNA_H
#define RNA_H

typedef double** partitionFunctionMatrix;

typedef double* partitionFunctionVector;


typedef struct {
  int length;
  double temperature;
  EnergyModel* energyModel;
  partitionFunctionVector Z;
  partitionFunctionMatrix Zb;
  partitionFunctionMatrix Z1;
  partitionFunctionMatrix Z2;
} RNA;




#endif
