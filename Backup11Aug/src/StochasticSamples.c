#include "StochasticSamples.h"
#include "EnergyFunctions.h"
#include "EnergyModel.h"
#include "RNA.h"
#include "AllowedPairs.h"
#include "PairIterator.h"
#include "PartitionFunction.h"
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

typedef int bool;


StochasticSamples* allocateStochasticSamples(RNA* strand, int numSamples) {

  StochasticSamples* newSamples = (StochasticSamples*) malloc(sizeof(StochasticSamples));
  newSamples->numSamples = numSamples;
  newSamples->structures = (Structure*) calloc(numSamples, sizeof(Structure));
  /* for(i = 0; i < numSamples; i++) { */
  /*   newSamples->structures[i] = (Structure) allocateStructure(strand); */
  /* }  */
  return newSamples;
}

Structure allocateStructure(RNA* strand) {
  return (Structure) calloc(strand->length, sizeof(int));
}

void freeStructure(Structure structure) {
  free(structure);
}

void pair(Structure s, int i, int j) {
  s[i] = j;
  s[j] = i;
}

void push(stackNode** s, int i, int j, int m) {
  struct stackNode* new_top;
  new_top = malloc(sizeof(struct stackNode));
  new_top->i = i;
  new_top->j = j;
  new_top->matrix = m;
  new_top->next = *s;
  *s = new_top;
}

void freeStochasticSamples(StochasticSamples* samples) {
  int i;
  for(i = 0; i < samples->numSamples; i++) {
    freeStructure(samples->structures[i]);
  } 
  free(samples->structures);
  free(samples);
}

void sample(RNA* strand, int samples) {
  if(strand->samples) {
    freeStochasticSamples(strand->samples);
  }

  strand->samples = allocateStochasticSamples(strand, samples);

  srand((unsigned int) strand->sequence);
  int numruns = 0, sample;
  for(sample =0; sample < samples; sample++) {
    while(!strand->samples->structures[sample]) { 
      numruns++;
      if(numruns > 100*samples) {
	exit(-1);
      }
      strand->samples->structures[sample] = sampleStructure(strand);
      if(!strand->samples->structures[sample]) { 
	//printf("Bad structure (sample: %d) (numruns: %d)\n", sample, numruns);
      }
    }
  }
}

Structure sampleStructure(RNA* strand) { 
  int i, j, k;
  double rnd, 
    multiA = strand->energyModel->multiA, 
    multiB = strand->energyModel->multiB, 
    multiC = strand->energyModel->multiC;
  bool found = FALSE;
  struct stackNode *stack, *top;
  double **Z = strand->partitionFunction->Z, 
    **Zb = strand->partitionFunction->Zb,
    **Z1 = strand->partitionFunction->Z1,
    **Z2 = strand->partitionFunction->Z2;
  double* scale = strand->energyModel->scale; 
  int len = strand->length;
  PairIterator* iteratorPlus1, *iterator, *iteratorMinus1;

  Structure sample = allocateStructure(strand);

  stack = NULL;
  j = len - 1;
  while (j > 0) {
    found = FALSE;
    rnd = (double) rand() / RAND_MAX * Z[0][j];
    if (rnd <= Z[0][j - 1] / scale[1]) {
      --j;
      continue;
    }
    else 
      rnd -= Z[0][j - 1] / scale[1];
    
    
    iterator = strand->allowedPairs->ji[j];
    for(k = start(iterator); hasNext(iterator); k = next(iterator)) {      
      if (rnd <= auPenalty(strand, k, j) * Zb[k][j] / scale[k]) { 
	push(&stack, k, j, 0);
	j = 0;
	found = TRUE;
	break;
      }
      else
	rnd -= auPenalty(strand, k, j) * Zb[k][j] / scale[k];
      
      if(k > 0) {
	if(rnd <= auPenalty(strand, k, j) * Z[0][k-1] * Zb[k][j]) {
	  push(&stack, k, j, 0);
	  j = k-1;
	  found = TRUE; 
	  break;
	} else 
	  rnd -= auPenalty(strand, k, j) * Z[0][k-1] * Zb[k][j];
	
	if(rnd <= auPenalty(strand, k, j) * Zb[k][j]/ scale[k-1] * ed5(strand, k, j)) {
	  //setDangle5(k, upst, dnst);
	  push(&stack, k, j, 0);
 	  j = k - 2;
	  found = TRUE;
	  break;
	} else 
	  rnd -= auPenalty(strand, k, j) * Zb[k][j] / scale[k-1] * ed5(strand, k, j);
      }
      
      if(k > 1) {
	if(rnd <= auPenalty(strand, k, j) * Z[0][k-2] * Zb[k][j] * ed5(strand, k, j)) {
	  //setDangle5(k, upst, dnst);
	  push(&stack, k, j, 0);
	  j = k-2;
	  found = TRUE;
	  break;
	} else 
	  rnd -= auPenalty(strand, k,j) * Z[0][k-2] * Zb[k][j] * ed5(strand, k, j);
      }
    } // j pairs for loop
    
    iteratorMinus1 = strand->allowedPairs->ji[j-1];    
    if(!found) {
      for(k = start(iteratorMinus1); hasNext(iteratorMinus1); k = next(iteratorMinus1))  {      
	
	if(rnd <= auPenalty(strand, k, j-1) * Zb[k][j-1] * ed3(strand, k, j-1) / scale[k]) {
	  //setDangle3(j-1, upst, dnst);
	  push(&stack, k, j-1, 0);
	  j = 0;
	  found = TRUE;
	  break;
	} else
	  rnd -= auPenalty(strand, k, j-1) * Zb[k][j-1] * ed3(strand, k, j-1) / scale[k];
	
	if(k > 0) {
	  if(rnd <= auPenalty(strand, k, j-1) * Z[0][k - 1] * Zb[k][j-1] * ed3(strand, k, j-1)) {
	    //setDangle3(j-1, upst, dnst);
	    push(&stack, k, j-1, 0);
	    j = k - 1;
	    found = TRUE;
	    break;
	  } else 
	    rnd -= auPenalty(strand, k, j-1) * Z[0][k - 1] * Zb[k][j-1] * ed3(strand, k, j-1);
	  
	  if(rnd <= auPenalty(strand, k, j-1) * Zb[k][j-1] * etstackm(strand, k, j-1) / scale[k-1]) {
	    //setDangle5(k, upst, dnst);
	    //setDangle3(j-1, upst, dnst);
	    push(&stack, k, j-1, 0);
	    j = 0;
	    found = TRUE;
	    break;
	  } else
	    rnd -= auPenalty(strand, k, j-1) * Zb[k][j-1] * etstackm(strand, k, j-1) / scale[k-1];
	}
	
	if(k > 1) {
	  if(rnd <= auPenalty(strand, k, j-1) * Z[0][k-2] * Zb[k][j-1] * etstackm(strand, k, j-2)) {
	    //setDangle5(k, upst, dnst);
	    //setDangle3(j-1, upst, dnst);
	    push(&stack, k, j-1, 0);
	    j = k - 2;
	    found = TRUE;
	  } else 
	    rnd -= auPenalty(strand, k, j-1) * Z[0][k-2] * Zb[k][j-1] * etstackm(strand, k, j-2);
	}
      } // j-1 pairs for loop
    }
    
    // in the case that there is nothing probably found treat as a draw for j -1
    if(!found) {
      free(sample);
      return 0;
    }
  }
  
  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      found = FALSE;
      if (top->matrix == 0) /* Zb */
	{
	  /* bp[i - 1] = j; */
	  /* bp[j - 1] = i; */
	  pair(sample, i, j);
	  rnd = (double) rand() / RAND_MAX * Zb[i][j];
	  if (rnd <= hairpinTerm(strand, i, j)) {
	    //hairpin, do nothing
	    found = TRUE; 
	  }
	  else if (rnd <= hairpinTerm(strand, i, j) + stackTerm(strand, i, j) * Zb[i + 1][j - 1])
	    { // stack
	      /* upst[i] = i; */
	      /* dnst[i - 1] = i + 1; */
	      /* upst[j - 1] = j - 1; */
	      /* dnst[j - 2] = j; */
	      push(&stack, i + 1, j - 1, 0);
	      found = TRUE;
	    }
	  else if (rnd <= hairpinTerm(strand, i, j) + stackTerm(strand, i, j) * Zb[i + 1][j - 1] + bulgeInternalTerm(strand, i, j)) {
	    // internal loop
	    push(&stack, i, j, 3);
	    found = TRUE;
	  }
	  else
	    {
	      // these are the multiloop terms
	      rnd -= hairpinTerm(strand, i, j) + stackTerm(strand, i, j) * Zb[i + 1][j - 1] + bulgeInternalTerm(strand, i, j);
	      rnd /= multiA * multiC * auPenalty(strand, i, j);
	      if(rnd <= Z2[i+1][j-1]) {
		push(&stack, i+1, j-1, 2);
		found = TRUE;
		continue;
	      } else
		rnd -= Z2[i+1][j-1];
	      
	      if (rnd <= ed5(strand, j, i)* multiB * Z2[i+1][j-2]) {
		//setDangle5(j, upst, dnst);
		push(&stack, i+1, j - 2, 2);
		found = TRUE;
		continue;
	      } else 
		rnd -= ed5(strand, j, i) * multiB * Z2[i+1][j-2];
	      
	      if(rnd <= ed3(strand, j, i) * multiB * Z2[i+2][j-1]) {
		//setDangle3(i, upst, dnst);
		push(&stack, i+2, j - 1, 2);
		found = TRUE;
		continue;
	      } else 
		rnd -= ed3(strand,j,i) * multiB * Z2[i+2][j-1];
	      
	      if(rnd <= etstackm(strand, j,i) * multiB * multiB *Z2[i+2][j-2]) {
		//setDangle5(j, upst, dnst);
		//setDangle3(i, upst, dnst);
		push(&stack, i+2, j-2, 2);
		found = TRUE;
		continue;
	      } else
		rnd -= etstackm(strand, j,i) * multiB * multiB * Z2[i+2][j-2];
	      
	      // should not get here
	      
	      if(!found) {
		//free(sample);
		//return 0;
		printf("Zb sampling not exhaustive, something is wrong\n");
		exit(-1);
	      }
	    }
	}
      else if (top->matrix == 1) /* Z1 */
	{
	  rnd = (double) rand() / RAND_MAX * Z1[i][j];
	  while (rnd <= multiB * Z1[i + 1][j] / scale[1])
	    {
	      i++;
	      rnd = (double) rand() / RAND_MAX * Z1[i][j];
	    }
	  rnd -= multiB * Z1[i + 1][j] / scale[1];
	  
	  iterator = strand->allowedPairs->ij[i];
	  for(k = start(iterator); hasNext(iterator); k = next(iterator)) { 
	    
	    if(rnd <= multiC * auPenalty(strand, i, k) * Zb[i][k] * scale[j-k]) {
	      push(&stack, i, k, 0);
	      found = TRUE;
	      break;
	    } else 
	      rnd -= multiC * auPenalty(strand, i, k) * Zb[i][k] * scale[j-k]; 
	    
	    if(k < j) {
	      if(rnd <= multiC * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+1][j]) {
		push(&stack, i, k, 0);
		push(&stack, k+1, j, 1);
		found = TRUE;
		break; 
	      } else 
		rnd -= multiC * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+1][j]; 
	    }
	    
	    
	    if(k < j) {
	      if(multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * scale[j-k - 1] * ed3(strand, i, k)){
		//setDangle3(k, upst, dnst);
		push(&stack, i, k, 0);
		found = TRUE;
		break; 
	      } else
		rnd -= multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * scale[j-k - 1] * ed3(strand, i, k);
	    }
	    
	    
	    if(k < j-1) {
	      if(rnd <= multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k)) {
		//setDangle3(k, upst, dnst);
		push(&stack, i, k, 0);
		push(&stack, k+2, j, 1);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	    }
	    
	  }
	  
	  iteratorPlus1 = strand->allowedPairs->ij[i + 1];
	  for(k = start(iteratorPlus1); hasNext(iteratorPlus1); k = next(iteratorPlus1)) {
	    if(k > j) break;
	    if(rnd <= multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * scale[j-k] * ed5(strand, i+1, k)){
	      //setDangle5(i+1, upst, dnst);
	      push(&stack, i+1, k, 0);
	      found = TRUE;
	      break;
	    } else 
	      rnd -= multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * scale[j-k] * ed5(strand, i+1, k);
	    
	    if(k < j) {
	      if (rnd <= multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k)) {
		//setDangle5(i+1, upst, dnst);
		push(&stack, i+1, k, 0);
		push(&stack, k+1, j, 1);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	    } 
	    
	    
	    if(k < j) {
	      if(rnd <= multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * scale[j-k - 1] * etstackm(strand, i+1, k)) {
		//setDangle5(i+1, upst, dnst);
		//setDangle3(k, upst, dnst);
		push(&stack, i+1, k, 0);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * scale[j-k - 1] * etstackm(strand, i+1, k);
	    }
	    
	    if(k < j - 1) {
	      if(rnd <= multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k)) { 
		//setDangle5(i+1, upst, dnst);
		//setDangle3(k, upst, dnst);
		push(&stack, i+1, k, 0);
		push(&stack, k+2, j, 1);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);
	    }
	  }
	  
	  if(!found) {
	    free(sample);
	    return 0;
	  }
	  
	} 
      else if (top->matrix == 2) /* Z2 */
	{
	  rnd = (double) rand() / RAND_MAX * Z2[i][j];
	  while (rnd <= multiB * Z2[i + 1][j] / scale[1])
	    {
	      i++;
	      rnd = (double) rand() / RAND_MAX * Z2[i][j];
	    }
	  rnd -= multiB * Z2[i + 1][j] / scale[1];
	  //Z2[i][j] = multiB * Z2[i+1][j] / scale;
	  
	  iterator = strand->allowedPairs->ij[i];	    
	  for(k = start(iterator); hasNext(iterator); k = next(iterator)) {
	    if(k > j) break;
	    if(k < j) {
	      // Z2[i][j] += multiC * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+1][j];	  
	      if(rnd <= multiC * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+1][j]) {
		push(&stack, i, k, 0);
		push(&stack, k+1, j, 1);
		found = TRUE;
		break; 
	      } else 
		rnd -= multiC * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+1][j]; 
	    }

	    if(k < j-1) {
	      // Z2[i][j] += multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	      if(rnd <= multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k)) {
		//setDangle3(k, upst, dnst);
		push(&stack, i, k, 0);
		push(&stack, k+2, j, 1);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * auPenalty(strand, i, k) * Zb[i][k] * Z1[k+2][j] * ed3(strand, i, k);
	    }
	  }
	  
	  iteratorPlus1 = strand->allowedPairs->ij[i+1];
	  for(k = start(iteratorPlus1); hasNext(iteratorPlus1); k = next(iteratorPlus1)) {	    	    	    
	    if(k > j) break;
	    if(k < j) {
	      // Z2[i][j] += multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	      if(rnd <= multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k)) {
		//setDangle5(i+1, upst, dnst);
		push(&stack, i+1, k, 0);
		push(&stack, k+1, j, 1);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+1][j] * ed5(strand, i+1, k);
	    }
	    if(k < j - 1) {
	      // Z2[i][j] += multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);
	      if(rnd <= multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k)) { 
		//setDangle5(i+1, upst, dnst);
		//setDangle3(k, upst, dnst);
		push(&stack, i+1, k, 0);
		push(&stack, k+2, j, 1);
		found = TRUE;
		break;
	      } else 
		rnd -= multiC * multiB * multiB * auPenalty(strand, i+1, k) * Zb[i+1][k] * Z1[k+2][j] * etstackm(strand, i+1, k);
	    }
	  }
	  
	  if(!found) {
	    free(sample);
	    return 0;
	  }
	  
	}
      else /* QBI */
	{
	  int d, ii, jj, done;
	  rnd = (double) rand() / RAND_MAX * bulgeInternalTerm(strand, i, j);
	  done = 0;
	  // for loop
	  for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - MAXLOOP && !done; --d)
	    for (ii = i + 1; ii < j - d && !done; ++ii)
	      {
		jj = d + ii;
		if (Zb[ii][jj] != 0)
		  {
		    if (rnd <= ebi(strand, i, j, ii, jj) * Zb[ii][jj])
		      {
			//setBI(i, j, ii, jj, upst, dnst);
			push(&stack, ii, jj, 0);
			done = 1;
			break;
		      }
		    else
		      rnd -= ebi(strand, i, j, ii, jj) * Zb[ii][jj];
		  }
	      }
	  if(!done) {
	    free(sample); 
	    return 0;
	  }
	}
      
    }
  return sample;
}
