#include "RNA.h"
#include "EnergyModel.h"
#include "EnergyFunctions.h"
#include "PartitionFunction.h"
#include <math.h>
#include <string.h>


// used to be g_misc[12]
#define HAIRPINBASE 37

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

double hairpinTerm(RNA* strand, int i, int j)
{
  double energy;
  int loopSize = j - i - 1;
  int k;
  int length = strand->length;
  EnergyModel* em = strand->energyModel;
  char* seq = strand->sequence;
  int* nSeq = strand->intSequence;

  if (loopSize < TURN)
    return 0.0;

  if (i <= length && length < j)
    return 0.0;
  else if (i > length)
    {
      i -= length;
      j -= length;
    }

  if (loopSize <= 30)
    energy = em->hairpinLoop[loopSize - 1];
  else
    energy = em->hairpinLoop[29] * pow(HAIRPINBASE, log((double) loopSize / 30)) / em->scale[loopSize - 30];

  
  if (loopSize > 3)
    energy *= em->tstack[nSeq[i]][nSeq[i+1]][nSeq[j-1]][nSeq[j]];
  else
    energy *= auPenalty(strand, i, j);
  
  // Markham did binary search here, not sure how much this will
  // effect performance, but requires sorting the data alphabetically
  // I believe.
  if (loopSize == 3)
    {
      //if ((loop = bsearch(g_seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
      for(k = 0; k < NUM_TRILOOPS; k++) {
	if(memcmp(em->triloop[k].loop, seq + i, 5)) {
	  energy *= em->triloop[k].energy;
	  break;
	}
      }
    }
  else if (loopSize == 4)
    {
      for(k = 0; k < NUM_TLOOPS; k++) {
	if(memcmp(em->tloop[k].loop, seq + i, 6)) {
	  energy *= em->tloop[k].energy;
	  break;
	}
      }
    }
  else if (loopSize == 6)
    {
      for(k = 0; k < NUM_HEXALOOPS; k++) {
	if(memcmp(em->hexaloop[k].loop, seq + i, 8)) {
	  energy *= em->hexaloop[k].energy;
	  break;
	}
      } 
    }

  /* /\* GGG *\/ */
  /* if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3) */
  /*   energy *= g_misc[8]; */

  /* /\* poly-C *\/ */
  /* if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1) */
  /*   energy *= g_misc[11]; */
  /* else */
  /*   { */
  /*     for (k = 1; k <= loopSize; ++k) */
  /* 	if (g_seq[i + k] != 1) */
  /* 	  return energy; */
  /*     energy *= pow(g_misc[9], loopSize) * g_misc[10]; */
  /*   } */

  return energy;
}

double bulgeInternalTerm(RNA* strand, int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;
  double** Zb = strand->partitionFunction->Zb;

  for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - MAXLOOP; --d)
    for (ii = i + 1; ii < j - d && ii <= strand->length; ++ii)
      {
	jj = d + ii;
	if(isCannonical(strand, ii, jj))
	  energy += ebi(strand, i, j, ii, jj) * Zb[ii][jj];
      }

  return energy;
}


double bulgeInternalTermPeriodic(RNA* strand, int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;
  double** Zb = strand->partitionFunction->Zb;
  for (d = j - i - 3; d >= 1 && d >= j - i - 2 - MAXLOOP; --d)
    for (ii = MAX(i + 1, strand->length + 1 - d); ii < j - d && ii <= strand->length; ++ii)
      {
	jj = d + ii;
	if (isCannonical(strand, ii, jj - strand->length))
	  energy += ebi(strand, jj - strand->length, ii, j - strand->length, i) * Zb[ii][jj];
      }

  return energy;
}


double auPenalty(RNA* strand, int i, int j) {
  int* nSeq = strand->intSequence;
  return strand->energyModel->auPenalty[nSeq[i]][nSeq[j]];
}

double stackTerm(RNA* strand, int i, int j)
{
  int length = strand->length;
  int* nSeq = strand->intSequence;
  if (i == length || j == length + 1)
    return 0.0;

  if (i > length)
    i -= length;
  if (j > length)
    j -= length;

  return strand->energyModel->stack[nSeq[i]][nSeq[i + 1]][nSeq[j - 1]][nSeq[j]];
}

double etstackm(RNA* strand, int i, int j)
{
  int* nSeq = strand->intSequence;

  return strand->energyModel->stack[nSeq[i]][nSeq[i + 1]][nSeq[j-1]][nSeq[j]];
}

double ed3(RNA* strand, int i, int j)
{
  int* nSeq = strand->intSequence;
  return strand->energyModel->dangle3[nSeq[j]][nSeq[i]][nSeq[j + 1]];
}

double ed5(RNA* strand, int i, int j)
{
  int* nSeq = strand->intSequence;
  return strand->energyModel->dangle5[nSeq[j]][nSeq[i]][nSeq[i - 1]];
}

double ebi(RNA* strand, int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  double loopEnergy;//, asPenalty;
  EnergyModel* em = strand->energyModel;
  int* nSeq = strand->intSequence;
  double scale = em->scale[1];

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > MAXLOOP)
    return 0.0;

  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return em->bulgeLoop[0] * em->stack[nSeq[i]][nSeq[j]][nSeq[ii]][nSeq[jj]] * scale * scale;
      else if (loopSize2 <= 30)
	return em->bulgeLoop[loopSize2 - 1] * auPenalty(strand, i, j) * auPenalty(strand, ii, jj);
      else
	return em->bulgeLoop[29] * pow(HAIRPINBASE, log((double) loopSize2 / 30)) / em->scale[loopSize2 - 30] * auPenalty(strand, i, j) * auPenalty(strand, ii, jj);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return em->bulgeLoop[0] * em->stack[nSeq[i]][nSeq[j]][nSeq[ii]][nSeq[jj]] * scale * scale;
      else if (loopSize1 <= 30)
	return em->bulgeLoop[loopSize1 - 1] * auPenalty(strand, i, j) * auPenalty(strand, ii, jj);
      else
	return em->bulgeLoop[29] * pow(HAIRPINBASE, log((double) loopSize1 / 30)) / em->scale[loopSize1 - 30] * auPenalty(strand, i, j) * auPenalty(strand, ii, jj);
    } /* Might leave these out (they are in the 2004 data tables, and we want 1999 for simplicity
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(nSeq[i], nSeq[j])][basePairIndex(nSeq[ii], nSeq[jj])][nSeq[i + 1]][nSeq[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(nSeq[i], nSeq[j])][basePairIndex(nSeq[ii], nSeq[jj])][nSeq[i + 1]][nSeq[j - 1]][nSeq[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(nSeq[jj], nSeq[ii])][basePairIndex(nSeq[j], nSeq[i])][nSeq[jj + 1]][nSeq[ii - 1]][nSeq[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(nSeq[i], nSeq[j])][basePairIndex(nSeq[ii], nSeq[jj])][nSeq[i + 1]][nSeq[j - 1]][nSeq[i + 2]][nSeq[j - 2]];
  else if ((loopSize1 == 2 && loopSize2 == 3) ||
	   (loopSize1 == 3 && loopSize2 == 2))
    return g_tstacki23[nSeq[i]][nSeq[j]][nSeq[i + 1]][nSeq[j - 1]] *
      g_tstacki23[nSeq[jj]][nSeq[ii]][nSeq[jj + 1]][nSeq[ii - 1]] / em->scale[7];
      */
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = em->internalLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = em->internalLoop[29] * pow(HAIRPINBASE, log((double) (loopSize1 + loopSize2) / 30)) / em->scale[loopSize1 + loopSize2 - 30];
      /* if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1)) */
      /* 	{ */
      /* 	  loopEnergy *= g_tstacki[nSeq[i]][nSeq[j]][0][0]; */
      /* 	  loopEnergy *= g_tstacki[nSeq[jj]][nSeq[ii]][0][0]; */
      /* 	} */
      /* else */
      /* 	{ */
      /* 	  loopEnergy *= g_tstacki[nSeq[i]][nSeq[j]][nSeq[i + 1]][nSeq[j - 1]]; */
      /* 	  loopEnergy *= g_tstacki[nSeq[jj]][nSeq[ii]][nSeq[jj + 1]][nSeq[ii - 1]]; */
      /* 	} */
      /* asPenalty = pow(g_misc[min3(4, loopSize1, loopSize2) - 1], abs(loopSize1 - loopSize2)); */
      /* if (asPenalty < g_misc[4]) */
      /* 	asPenalty = g_misc[4]; */
      /* loopEnergy *= asPenalty; */

      return loopEnergy;
    }
}
