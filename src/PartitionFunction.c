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
  ret->filledP = 0;
  ret->filledZbZ1Z2extended = 0;
  ret->Z = (double**) malloc(2 * length * sizeof(double*));
  ret->Zb = (double**) malloc(2 * length * sizeof(double*));
  ret->Z1 = (double**) malloc(2 * length * sizeof(double*));
  ret->Z2 = (double**) malloc(2 * length * sizeof(double*));
  ret->P = (double**) malloc(length * sizeof(double*));
  int i;
  for(i=0; i < 2*length; i++) {
    ret->Z[i] = (double*) calloc(2 * length,  sizeof(double));
    ret->Zb[i] = (double*) calloc(2 * length, sizeof(double));
    ret->Z1[i] = (double*) calloc(2 * length, sizeof(double));
    ret->Z2[i] = (double*) calloc(2 * length, sizeof(double));
    if(i < length) 
      ret->P[i] = (double*) calloc(length, sizeof(double));
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

// this just fills the first row and last column
void fillZ(RNA* strand) {
  int i,j, k, len = strand->length;
  double** Z = strand->partitionFunction->Z;
  double** Zb = strand->partitionFunction->Zb;
  double* scale = strand->energyModel->scale;
  PairIterator* iterator;
  PairIterator* iteratorMinus1;

  Z[0][0] = 0.0;

  for (j = 1; j < len; ++j) {

    Z[0][j] = Z[0][j - 1] / scale[1];
    iterator = strand->allowedPairs->ji[j];
    for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
      if(k < 0) break;
      Z[0][j] += auPenalty(k, j) * Zb[k][j] / scale[k - 1];
      if(k > 1) {
	Z[0][j] += auPenalty(k, j) * Z[0][k-1] * Zb[k][j];
	Z[0][j] += auPenalty(k, j) * Zb[k][j] / scale[k-2] * ed5(strand, k, j);
	if(k > 2) {
	  Z[0][j] += auPenalty(k, j) * Z[0][k-2] * Zb[k][j] * ed5(strand, k, j);
	}
      }
    }
      
    iteratorMinus1 = strand->allowedPairs->ji[j-1];
    for(k = start(iteratorMinus1); hasNext(iteratorMinus1); k = next(iteratorMinus1)) {
      if(k < 0) break;
      Z[0][j] += auPenalty(k, j-1) * Zb[k][j-1] * ed3(strand, k, j-1) / scale[k-1];
      if(k > 1) {
	Z[0][j] += auPenalty(k, j-1) * Z[0][k-1] * Zb[k][j-1] * ed3(strand, k, j-1);
	Z[0][j] += auPenalty(k, j-1) * Zb[k][j-1] * etstackm(strand, k, j-1) / scale[k-2];
	if(k > 2) {
	  Z[0][j] += auPenalty(k, j-1) * Z[0][k-2] * Zb[k][j-1] * etstackm(strand, k, j-1);
	}
      }
    
    }
  }


  for(i = len - 2; i >= 0; i--) {
    Z[i][len-1] = Z[i + 1][len-1] / scale[1];
    
    iterator = strand->allowedPairs->ij[i];
    for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
      Z[i][len-1] += auPenalty(i, k) * Zb[i][k]/scale[len-1 - k];
      if(k < len-1) {
	Z[i][len-1] += auPenalty(i, k) * Zb[i][k] * Z[k+1][len-1]; 
	Z[i][len-1] += auPenalty(i, k) * Zb[i][k] / scale[len-1 - k - 1] * ed3(strand,i,k);
	if(k < len-1-1) {
	  Z[i][len-1] += auPenalty(i, k) * Zb[i][k] * Z[k+2][len-1] * ed3(strand, i, k);
	}
      }
    }
    
    iteratorMinus1 = strand->allowedPairs->ij[i+1];
    for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
      Z[i][len-1] += auPenalty(i+1, k) * Zb[i+1][k] / scale[len-1 - k] * ed5(strand, i+1, k);
      if(k < len-1) {
	Z[i][len-1] += auPenalty(i+1, k) * Zb[i+1][k] * Z[k+1][len-1] * ed5(strand, i+1, k); 
	Z[i][len-1] += auPenalty(i+1, k) * Zb[i+1][k] / scale[len-1 - k - 1] * etstackm(strand, i + 1, k);
	if(k < len-1-1) {
	  Z[i][len-1] += auPenalty(i+1, k) * Zb[i+1][k] * Z[k+2][len-1] * etstackm(strand, i + 1, k);
	}
      }
    }
  }

  strand->partitionFunction->filledZ = 1;
}


void fillExtendedZbZ1Z2(RNA* strand) {
  int i,j,k,len = strand->length;

  double** Z = strand->partitionFunction->Z;
  double** Zb = strand->partitionFunction->Zb;
  double** Z2 = strand->partitionFunction->Z2;
  double** Z1 = strand->partitionFunction->Z1;


  double multiA = 1;//strand.energyModel.multiA;
  double multiB = 1;//strand.energyModel.multiB;
  double multiC = 1;//strand.energyModel.multiC;
  double* scale = strand->energyModel->scale;
  double* bscale = strand->energyModel->bscale;
  PairIterator* iterator;
  PairIterator* iteratorPlus1;

  for(j = len + 1; j < 2 * len; j++) {
    for(i = len - 2; i > j - len; i--) {
      double au = auPenalty(i, j - len);
      // no hairpin term, not allowed in crossing 
      Zb[i][j] = stackTerm(strand, i, j-len) * Zb[i+1][j-1] + bulgeInternalTerm(strand, i, j - len);

      Zb[i][j] += multiA * multiC * au * (Z2[i+1][j-1]
      					  + ed5(strand, j - len, i) * multiB *Z2[i+1][j-2]
      					  + ed3(strand, j - len,i) * multiB * Z2[i+2][j-1]
      					  + etstackm(strand, j - len,i) * multiB * multiB * Z2[i+2][j-2]);


      Zb[i][j] += au * (Z[i+1][len - 1] + 1/scale[len - 1 - i]) *
	(Z[0][j - 1 - len] + 1/scale[j-len])/scale[2];
      
      if(j > len + 1)  {
	Zb[i][j] += au * ed5(strand, j - len, i) * (Z[i + 1][len - 1] + 1/scale[len - 1 - i]) * 
	  (Z[0][j - 2 - len] + 1/scale[j-len-1])/scale[2]; 
      }
      if(i < len - 2) {
	Zb[i][j] += au * ed3(strand, j - len, i) * (Zb[i + 2][len - 1] + 1/scale[len - 2 - i]) * 
	  (Zb[0][j - 1 - len] + 1/scale[j - len])/scale[2];
      }
      if(i < len - 2 && j > len + 1 ) {
	Zb[i][j] += au * etstackm(strand, j - len, i) * (Zb[i+2][len-1] + 1/scale[len - 2 - i]) * 
	  (Zb[0][j - 2 - len] + 1/scale[j-len-1])/scale[2];
      }


      Z2[i][j] = multiB * Z2[i+1][j] / scale[0];
      Z1[i][j] = multiB * Z1[i+1][j] / scale[0];	

      // k is on the same side as i *********************************
      // in this case Z1 is the same as Z2 because no completely empty
      // sections are allowed between k and j since they are on
      // different sides, this changes in the next section
      iterator = strand->allowedPairs->ij[i];
      for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
	if(k < len - 1) {
	  Z1[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * Z1[k+1][j]; 
	  Z2[i][j] += multiC * auPenalty(i, k) * Zb[i][k] * Z1[k+1][j];	  
	  if(k < len - 2) {
	    Z1[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	    Z2[i][j] += multiC * multiB * auPenalty(i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	  }
	}
      }

      iteratorPlus1 = strand->allowedPairs->ij[i+1];
      for(k = start(iteratorPlus1); hasNext(iteratorPlus1); k = next(iteratorPlus1)) {
	if(k < len - 1) { 
	  Z1[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	  Z2[i][j] += multiC * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	  if(k < len + 2) {
	    Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);	  
	    Z2[i][j] += multiC * multiB * multiB * auPenalty(i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);
	  }
	}
      }

      // k is on the same side as j ********************************
      iterator = strand->allowedPairs->ji[i];
      for(k = start(iterator) + len; hasNext(iterator); k = next(iterator) + len) {
	if(k > j) continue;
	Z1[i][j] += multiC * auPenalty(i, k - len) * Zb[i][k] * bscale[j-k];	  
	if(k < j) { 
	  Z1[i][j] += multiC * auPenalty(i, k - len) * Zb[i][k] * Z1[k+1 - len][j - len]; 
	  Z1[i][j] += multiC * multiB * auPenalty(i, k - len) * Zb[i][k] * bscale[j-(k+1)] * ed3(strand, i, k - len);
	  Z2[i][j] += multiC * auPenalty(i, k - len) * Zb[i][k] * Z1[k+1 - len][j - len];	  
	  if(k < j - 1) {
	    Z1[i][j] += multiC * multiB * auPenalty(i, k - len) * Zb[i][k] * Z1[k+2 - len][j - len] * ed3(strand, i, k - len);
	    Z2[i][j] += multiC * multiB * auPenalty(i, k - len) * Zb[i][k] * Z1[k+2 - len][j - len] * ed3(strand, i, k - len);
	  }
	}
      }

      iteratorPlus1 = strand->allowedPairs->ji[i+1];
      for(k = start(iteratorPlus1) + len; hasNext(iteratorPlus1); k = next(iteratorPlus1) + len) {
	if(k > j) continue;
	Z1[i][j] += multiC * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * bscale[j-k] * ed5(strand, i+1, k - len);
	if(k < j) { 
	  Z1[i][j] += multiC * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+1 - len][j - len] * ed5(strand, i+1, k - len);
	  Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * bscale[j- (k+1)] * etstackm(strand, i+1, k - len);
	  Z2[i][j] += multiC * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+1 - len][j - len] * ed5(strand, i+1, k - len);
	  if(k < j - 1) {
	    Z1[i][j] += multiC * multiB * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+2 - len][j - len] * etstackm(strand, i+1, k - len);	  
	    Z2[i][j] += multiC * multiB * multiB * auPenalty(i+1, k - len) * Zb[i+1][k] * Z1[k+2 - len][j - len] * etstackm(strand, i+1, k - len);
	  }
	}
      }  
    } 
  }

  strand->partitionFunction->filledZbZ1Z2extended = 1;
}


void fillP(RNA* strand) {
  int i, j, len = strand->length;
  double** P = strand->partitionFunction->P;
  double** Zb = strand->partitionFunction->Zb;
  double** Z = strand->partitionFunction->Z;
  double* scale = strand->energyModel->scale;
  for(i =  0; i < len; i++) {
    for(j = i + TURN + 1; j < len; j++) {
      if(Zb[i][j] == 0)
	continue;
      P[i][j] = Zb[i][j] * Zb[j][i + len] * scale[2] / Z[0][len - 1];
      if(P[i][j] > 1.001 || P[i][j] < -.001) {
	printf("Warning: P(%d,%d) = %g\n", i,j,P[i][j]);
      }
    }
  }
  strand->partitionFunction->filledP = 1;
}
