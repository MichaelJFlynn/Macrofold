// Written by Mike Flynn, 2015-2016

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

	PairIterator* iterator = strand->allowedPairs->ij[i];
	PairIterator* iteratorPlus1 = strand->allowedPairs->ij[i+1];


	Zb[i][j] = hairpinTerm(strand, i, j) + stackTerm(strand, i,j) * Zb[i + 1][j - 1] + bulgeInternalTerm(strand, i, j);	
	Zb[i][j] += multiA * multiC * au * (Z2[i+1][j-1] 
					  + ed5(strand, j, i) * multiB *Z2[i+1][j-2]
					  + ed3(strand, j,i) * multiB * Z2[i+2][j-1] 
					  + etstackm(strand, j,i) * multiB * multiB * Z2[i+2][j-2]);
	

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
}


void fillExtendedZbZ1Z2(RNA* strand) {
  int i,j,k,len = strand->length;

  double** Zb = strand->partitionFunction->Zb;
  double** Z2 = strand->partitionFunction->Z2;
  double** Z1 = strand->partitionFunction->Z1;


  double multiA = 1;//strand.energyModel.multiA;
  double multiB = 1;//strand.energyModel.multiB;
  double multiC = 1;//strand.energyModel.multiC;
  double scale = strand->energyModel->scale[1];
  double* bscale = strand->energyModel->bscale;
  PairIterator* iterator;
  PairIterator* iteratorPlus1;

  for(j = len + 1; j < 2 * len; j++) {
    for(i = len - 2; i > j - len; i++) {
      double au = auPenalty(i, j - len);

      // no hairpin term, not allowed in crossing 
      Zb[i][j] = stackTerm(strand, i, j) * Zb[i+1][j-1] + bulgeInternalTerm(strand, i, j);
      Zb[i][j] += multiA * multiC * au * (Z2[i+1][j-1] 
					  + ed5(strand, j, i) * multiB *Z2[i+1][j-2]
					  + ed3(strand, j,i) * multiB * Z2[i+2][j-1]
					  + etstackm(strand, j,i) * multiB * multiB * Z2[i+2][j-2]);
      
      Zb[i][j] += au * (Z[i+1][len - 1] + 1/strand->scale[len - 1 - i]) *
	(Z[0][j - 1 - len] + 1/strand->scale[j-len])/strand->scale[2];
      
      if(j > len + 1)  {
	Zb[i][j] += au * ed5(strand, j - len, i) * (Z[i + 1][len - 1] + 1/strand->scale[len - 1 - i]) * 
	  (Z[0][j - 2 - len] + 1/strand->scale[j-len-1])/strand->scale[2]; 
      }
      if(i < len - 2) {
	Zb[i][j] += au * ed3(strand, j - len, i) * (Zb[i + 2][len - 1] + 1/strand->scale[len - 2 - i]) * 
	  (Zb[0][j - 1 - len] + 1/strand->scale[j - len])/strand->scale[2];
      }
      if(i < len - 2 && j > len + 1 ) {
	Zb[i][j] += au * etstackm(strand, j - len, i) * (Zb[i+2][len-1] + 1/strand->scale[len - 2 - i]) * 
	  (Zb[0][j - 2 - len] + 1/strand->scale[j-len-1])/strand->scale[2];
      }

      Z2[i][j] = multiB * Z2[i+1][j] / scale;
      Z1[i][j] = multiB * Z1[i+1][j] / scale;	

      // k is on the same side as i *********************************
      // in this case Z1 is the same as Z2 because no completely empty
      // sections are allowed between k and j since they are on
      // different sides, this changes in the next section
      iterator = strand->allowedPairs->ij[i];
      for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
	Z1[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * Z1[k+1][j]; 
	Z2[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * Z1[k+1][j];	  
	Z1[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	Z2[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
      }

      iteratorPlus1 = strand->allowedPairs->ij[i+1];
      for(k = start(iteratorPlus1); hasNext(iteratorPlus1); k = next(iteratorPlus1)) {
	Z1[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	Z2[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);	  
	Z2[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);
      }
      
      // k is on the same side as j ********************************
      iterator = strand->allowedPairs->ji[i];
      for(k = start(iterator) + len; hasNext(iterator); k = next(iterator) + len) {
	if(k > j) continue;
	Z1[i][j] += multiC * auPenalty(i, k - len) * Zb[i][k] * bscale[j-k];	  
	if(k < j) { 
	  Z1[i][j] += multiC * auPenalty(i, k - len) * Zb[i][k] * Z1[k+1][j]; 
	  Z1[i][j] += multiC * multiB * auPenalty(i, k - len) * Zb[i][k] * bscale[j-(k+1)] * ed3(strand, i, k - len);
	  Z2[i][j] += multiC * auPenalty(i, k - len) * Zb[i][k] * Z1[k+1][j];	  
	  if(k < j - 1) {
	    Z1[i][j] += multiC * multiB * auPenalty(i, k - len) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k - len);
	    Z2[i][j] += multiC * multiB * auPenalty(i, k - len) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k - len);
	  }
	}
      }
      
      iteratorPlus1 = strand->allowedPairs->ji[i+1];
      for(k = start(iteratorPlus1) + len; hasNext(iteratorPlus1); k = next(iteratorPlus1) + len) {
	if(k > j) continue;
	Z1[i][j] += multiC * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * bscale[j-k] * ed5(strand, i+1, k - len);
	if(k < j) { 
	  Z1[i][j] += multiC * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k - len);
	  Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * bscale[j- (k+1)] * etstackm(strand, i+1, k - len);
	  Z2[i][j] += multiC * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k - len);
	  if(k < j - 1) {
	    Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k - len);	  
	    Z2[i][j] += multiC * multiB * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k - len);
	  }
	}
      }            
    } 
  }
}


void fillP(RNA* strand) {
  int i, j;
  for(i =  0; i < g_len; i++) {
    for(j = i + TURN + 1; j < g_len; j++) {
      if(Zb[i][j] == 0)
	continue;
      P(i,j) = Zb[i,j] * Zb[j][i + g_len] * g_scalen[2] / Z[g_len - 1];
      P -= Q0(i,j) * hairpinTerm(strand, i , j) / Z[g_len - 1];
      if(P(i,j) > 1.001 || P(i,j) < -.001) {
	printf("Warning: P(%d,%d) = %g\n", i,j,P(i,j));
      }
    }
  }
}
