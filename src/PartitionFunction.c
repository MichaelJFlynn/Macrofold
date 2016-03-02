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

  for(j = len; j < 2 * len; j++) {
    for(i = len - 1; i > j - len; i++) {
      double au = auPenalty(i, j - len);
      if(i < len - 1 && j > len) {
	Zb[i][j] = stackTerm(strand, i, j) * Zb[i+1][j-1] + bulgeInternalTerm(strand, i, j);
	Zb[i][j] += multiA * multiC * au * (Z2[i+1][j-1] 
					  + ed5(strand, j, i) * multiB *Z2[i+1][j-2]
					  + ed3(strand, j,i) * multiB * Z2[i+2][j-1] 
					  + etstackm(strand, j,i) * multiB * multiB * Z2[i+2][j-2]);
	/* for(k = i + TURN + 3; k < len; k++) { */
	/*   Zb[i][j] = multiA * multiC * au * Q(i+1, k-1) * Q1(k, j-1); */
	/*   Zb[i][j] = multiA * multiB * multiC * au * ed5(strand, j - len, i) * Q(i + 1, k - 1) * Q1(k, j - 2); */
	/*   Zb[i][j] = multiA * multiB * multiC * au * ed3(strand, j - len, i) * Q(i + 2, k - 1) * Q!(k, j - 1); */
	/*   Zb[i][j] = multiA * multiB * multiB * multiC * etstackm(strand, j - len, i)  Q(i+2, k -1) * Q1(k, j - 2); */
	/* } */
	/* for(k = len + 1; k < j - TURN - 1; k++) { */
	/*   Zb[i][j] = multiA * multiC * au * Q(i+1, k-1) * Q1(k - g_len, j - 1 - g_len); */
	/*   Zb[i][j] = multiA * multiB * multiC * au * ed5(strand, j - len, i) * Q(i+1, k-1) * Q1(k - g_len, j - 2 - g_len); */
	/*   if(i  < len - 1) { */
	/*     Zb[i][j] = multiA * multiB * multiC * au * ed3(strand, j - len, i) * Q(i+2, k - 1) * Q1(k - g_len, j - 1 - g_len); */
	/*     Zb[i][j] = multiA * multiB * multiB * multiC * au * etstackm(strand, j-g_len, i) * Q(i + 2, k-1) * Q1(k - len, j - 2 - g_len); */
	/*   } */
	/* } */
      }
      else
	Zb[i][j] = 0.0;

      Zb[i][j] += au * (Z[j-g_len - 1] + 1 / g_scalen[j-g_len -1]) / g_scalen[g_len - i] / g_scalen[2];
      if(j > g_len + 1)  {
	Zb[i][j] += au * ed5(strand, j - g_len, i) * (Q3(i + 1) + 1/g_scalen[j - g_len - 2]) / g_scalen[2];
      }
      if(i < g_len) {
	Zb[i][j] += au * ed3(strand, j - g_len, i) * (Q3(i + 2) + 1/g_scalen[j-g_len-1]) / g_scalen[2];
      }
      if( j > g_len + 1 && i < g_len) {
	Zb[i][j] += au * etstackm(strand, j - g_len, i) * (Q3(i + 2) + 1/g_scalen[j-g_len - 2]) / g_scalen[2];
      }


      Q1(i,j) = g_multi[2] * au * Zb[i][j];
      if(i < g_len) {
	Q1(i, j) += g_multi[1] * g_multi[2] * auPenalty(i + 1, j - g_len) * ed5(strand, i + 1, j - g_len) * Zb[i + 1][j];
      }
      if(j < g_len + 1) {
	Q1(i, j) += g_multi[1] * g_multi[2] * auPenalty(i, j - 1 - g_len) * ed5(strand,i, j - 1 - g_len) * Zb[i][j - 1];
      }
      if( i < g_len && j > g_len + 1) {
	Q1(i,j) += g_multi[1] * g_multi[2] * auPenalty(i + 1, j - 1 - g_len) * etstackm(strand, i + 1, j - 1 - g_len) * Zb[i+1][j-1];
      } 
      if(j > g_len + 1) {
	Q1(i,j) += g_multi[1] * Q1(i, j - 1)/ g_scale;
      }

      Q(i,j) = Q1(i,j);
      if(i < g_len) {
	Q(i, j) += g_multi[1] * Q1(i + 1, j) / g_scale;
      }
      for(k = i+2, k <= g_len; k++) {
	Q(i,j) += (Q(i,k-1) + g_bscalen[k-i]) * Q1(k, j);
      }
      for(k = g_len + 2; k < j - TURN; k++) {
	Q(i, j) += Q(i, k - 1) * Q1(k - g_len, j - g_len);
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
