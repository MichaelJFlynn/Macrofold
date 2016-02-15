// Written by Mike Flynn, 2015

#include "EnergyFunctions.h"
#include "RNA.h"
#include "EnergyModel.h"
#include "AllowedPairs.h"
#include "PairIterator.h"
#include "PartitionFunction.h"

#define TURN 3
#define auPenalty(i, j) 1

PartitionFunction* allocatePartitionFunction(int length) {
  PartitionFunction* ret = (PartitionFunction*) malloc(sizeof(PartitionFunction));
  ret->filledZbZ1Z2 = 0;
  ret->filledZ = 0;
  ret->Z = (double*) malloc(length * sizeof(double));
  ret->Zb = (double**) malloc(length * sizeof(double*));
  ret->Z1 = (double**) malloc(length * sizeof(double*));
  ret->Z2 = (double**) malloc(length * sizeof(double*));
  
  int i, j;
  for(i=0; i < length; i++) {
    ret->Zb[i] = (double*) malloc(length * sizeof(double));
    ret->Z1[i] = (double*) malloc(length * sizeof(double));
    ret->Z2[i] = (double*) malloc(length * sizeof(double));
    ret->Z[i] = 0;
    for(j=0; j<length; j++) {
      ret->Zb[i][j] = 0;
      ret->Z1[i][j] = 0;
      ret->Z2[i][j] = 0;
    }
  }

  return ret;
}

void freePartitionFunction(RNA* strand) {
  int i;  
  PartitionFunction* ret = strand->partitionFunction;
  for(i = 0; i < strand->length; i++) {
    free(ret->Zb[i]);
    free(ret->Z1[i]);
    free(ret->Z2[i]);
  }
  free(ret->Z);
  free(ret->Zb);
  free(ret->Z1);
  free(ret->Z2);
  free(ret);
}

void fillZbZ1Z2(RNA* strand) {
  int i, j, k;  

  //  double* Z = strand->pfData.Z;
  double** Zb = strand->partitionFunction->Zb;
  double** Z2 = strand->partitionFunction->Z2;
  double** Z1 = strand->partitionFunction->Z1;


  double multiA = 1;//strand.energyModel.multiA;
  double multiB = 1;//strand.energyModel.multiB;
  double multiC = 1;//strand.energyModel.multiC;
  double scale = strand->energyModel->scale[1];
  double* bscale = strand->energyModel->bscale;

  /*
    Recurrence relations:
    Z(i, i) = 1
    Z(i, j) = \sum_k Z^b(i, k)Z(i+1, k) + Z(i+1, j)
    
    Z^b(i, i) = 0
    Z^b(i, j) = Zh(i, j) + Zs(i,j)Z^b(i+1, j-1) + QBI(i,j)
    + e^{-\beta a} Z^2(i+1, j-1) -- dangles?

    Z2(i, i) = 0
    Z2(i, j) = e^{-\beta c} \sum_k Z^b(i, k) Z1(k+1, j) 
    + e^{-\beta b} Z2(i+1, j)

    Z1(i, i) = 0
    Z1(i, j) = e^{-\beta c} \sum_k Z^b(i, k)[ Z1(k+1, j) + Zempty(k+1, j)]
    + e^{-\beta b} Z1(i+1, j)
   */

  for (j = TURN + 1; j < strand->length; ++j)
    for (i = j - TURN - 1; i >= 0; --i)
      {
	double au = auPenalty(i, j);
	//double Qmulti = 0;
	double Zmulti = 0; 
	PairIterator* iterator = strand->allowedPairs->ij[i];
	PairIterator* iteratorPlus1 = strand->allowedPairs->ij[i+1];

	if(Zb[i][j] != 0.0) {
	  Zb[i][j] = hairpinTerm(strand, i, j) + stackTerm(strand, i,j) * Zb[i + 1][j - 1] + bulgeInternalTerm(strand, i, j);	
	  Zmulti += multiA * multiC * au * (Z2[i+1][j-1] 
					    + ed5(strand, j, i) * multiB *Z2[i+1][j-2]
					    + ed3(strand, j,i) * multiB * Z2[i+2][j-1] 
					    + etstackm(strand, j,i) * multiB * multiB * Z2[i+2][j-2]);
	  Zb[i][j] += Zmulti;
	}

	Z2[i][j] = multiB * Z2[i+1][j] / scale;
	Z1[i][j] = multiB * Z1[i+1][j] / scale;	
	for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
	  if(k > j) break;
	  Z1[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * bscale[j-k];	  
	  if(k < j) { 
	    Z1[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * Z1[k+1][j]; 
	    Z1[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * bscale[j-(k+1)] * ed3(strand, i, k);
	    Z2[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * Z1[k+1][j];	  
	    if(k < j - 1) {
	      Z1[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	      Z2[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	    }
	  }
	}

	for(k = start(iteratorPlus1); hasNext(iteratorPlus1); k = next(iteratorPlus1)) {
	  if(k > j) break;
	  Z1[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * bscale[j-k] * ed5(strand, i+1, k);
	  if(k < j) { 
	    Z1[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	    Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * bscale[j- (k+1)] * etstackm(strand, i+1, k);
	    Z2[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	    if(k < j - 1) {
	      Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);	  
	      Z2[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);
	    }
	  }
	}
	
      }

  strand->partitionFunction->filledZbZ1Z2 = 1;
}


void fillZ(RNA* strand) {
  int j, k;
  double* Z = strand->partitionFunction->Z;
  double** Zb = strand->partitionFunction->Zb;
  double* scale = strand->energyModel->scale;
  PairIterator* iterator;
  PairIterator* iteratorMinus1;

  Z[0] = 0.0;
  Z[strand->length - 1] = 0.0;

  for (j = 1; j < strand->length; ++j) {

    Z[j] = Z[j - 1] / scale[1];
    iterator = strand->allowedPairs->ji[j];
    for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
      if(k < 0) break;
      Z[j] += auPenalty(k, j) * Zb[k][j] / scale[k - 1];
      if(k > 1) {
	Z[j] += auPenalty(k, j) * Z[k-1] * Zb[k][j];
	Z[j] += auPenalty(k, j) * Zb[k][j] / scale[k-2] * ed5(strand, k, j);
	if(k > 2) {
	  Z[j] += auPenalty(k, j) * Z[k-2] * Zb[k][j] * ed5(strand, k, j);
	}
      }
    }
      
    iteratorMinus1 = strand->allowedPairs->ji[j-1];
    for(k = start(iteratorMinus1); hasNext(iteratorMinus1); k = next(iteratorMinus1)) {
      if(k < 0) break;
      Z[j] += auPenalty(k, j-1) * Zb[k][j-1] * ed3(strand, k, j-1) / scale[k-1];
      if(k > 1) {
	Z[j] += auPenalty(k, j-1) * Z[k-1] * Zb[k][j-1] * ed3(strand, k, j-1);
	Z[j] += auPenalty(k, j-1) * Zb[k][j-1] * etstackm(strand, k, j-1) / scale[k-2];
	if(k > 2) {
	  Z[j] += auPenalty(k, j-1) * Z[k-2] * Zb[k][j-1] * etstackm(strand, k, j-1);
	}
      }
    
    }
  }
  strand->partitionFunction->filledZ = 1;
  /* for (i = strand.length - 1; i >= 1; --i) { */
  /*   Z(i, strand.length) = Z(i + 1, strand.length) / scale; */

  /*   for(index=0; pfc.ij[i][index] != 0 && pfc.ij[i][index] <= strand.length; ++index) { */
  /*     k = pfc.ij[i][index]; */
  /*     Z(i,strand.length) += auPenalty(i, k) * Zb(i, k)/scalen[strand.length - k]; */
  /*     if(k < strand.length) { */
  /* 	Z(i,strand.length) += auPenalty(i, k) * Zb(i, k) * Z(k+1, strand.length); */
  /* 	Z(i, strand.length) += auPenalty(i, k) * Zb(i, k) / scalen[strand.length - k - 1] * ed3(i,k); */
  /* 	if(k < strand.length-1) { */
  /* 	  Z(i,strand.length) += auPenalty(i, k) * Zb(i, k) * Z(k+2, strand.length) * ed3(i, k); */
  /* 	} */
  /*     } */
  /*   } */

  /*   for(index=0; pfc.ij[i+1][index] != 0 && pfc.ij[i+1][index] <= strand.length; ++index) { */
  /*     k = pfc.ij[i+1][index]; */
  /*     Z(i, strand.length) += auPenalty(i+1, k) * Zb(i+1, k) / scalen[strand.length - k] * ed5(i+1, k); */
  /*     if(k < strand.length) { */
  /* 	Z(i,strand.length) += auPenalty(i+1, k) * Zb(i+1, k) * Z(k+1, strand.length) * ed5(i+1, k); */
  /* 	Z(i,strand.length) += auPenalty(i+1, k) * Zb(i+1, k) / scalen[strand.length - k - 1] * etstackm(i + 1, k); */
  /* 	if(k < strand.length-1) { */
  /* 	  Z(i,strand.length) += auPenalty(i+1, k) * Zb(i+1, k) * Z(k+2, strand.length) * etstackm(i + 1, k); */
  /* 	} */
  /*     } */
  /*   } */
  /* } */
}
