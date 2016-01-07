
double HairpinTerm(int i, int j)
{
  double energy;
  int loopSize = j - i - 1;
  int k;

  if (loopSize < TURN)
    return 0.0;

  if (i <= g_len && g_len < j)
    return 0.0;
  else if (i > g_len)
    {
      i -= g_len;
      j -= g_len;
    }

  if (loopSize <= 30)
    energy = g_hairpinLoop[loopSize - 1];
  else
    energy = g_hairpinLoop[29] * pow(g_misc[12], log((double) loopSize / 30)) / g_scalen[loopSize - 30];

  if (loopSize > 3)
    energy *= g_tstackh[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
  else
    energy *= auPenalty(i, j);

  if (loopSize == 3)
    {
      struct triloop* loop;
      if (numTriloops)
	if ((loop = bsearch(g_seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  energy *= loop->energy;
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = bsearch(g_seq + i, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  energy *= loop->energy;
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = bsearch(g_seq + i, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  energy *= loop->energy;
    }

  /* GGG */
  if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3)
    energy *= g_misc[8];

  /* poly-C */
  if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1)
    energy *= g_misc[11];
  else
    {
      for (k = 1; k <= loopSize; ++k)
	if (g_seq[i + k] != 1)
	  return energy;
      energy *= pow(g_misc[9], loopSize) * g_misc[10];
    }

  return energy;

}
/* TODO:
 * g_stack table [check] 
 * g_scalen table [check]
 * g_tstack table
 * g_tstackh
 * g_tstackm
 * g_tstacki
 * g_dangle5 table
 * g_dangle3 table
 * g_misc table
 * g_hairpinLoop table [check]
 * g_bulgeLoop [check]
 * g_internalLoop [check]
 * g_sint2
 * g_asint1x2
 * g_sint4
 * g_tstacki23
 * g_tloop [check]
 * g_triloop [check]
 * g_hexaloop [check]
 */

double BulgeInternalTerm(int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;

  for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop; --d)
    for (ii = i + 1; ii < j - d && ii <= g_len; ++ii)
      {
	jj = d + ii;
	if (Zb(ii, jj) != 0.0)
	  energy += Ebi(i, j, ii, jj) * Zb(ii, jj);
      }

  return energy;
}


double StackTerm(int i, int j)
{
  if (i == g_len || j == g_len + 1)
    return 0.0;

  if (i > g_len)
    i -= g_len;
  if (j > g_len)
    j -= g_len;

  return g_stack[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
}

double Etstackm(int i, int j)
{
  return g_tstackm[g_seq[j]][g_seq[i]][g_seq[j + 1]][g_seq[i - 1]];
}

double Ed3(int i, int j)
{
  return g_dangle3[g_seq[j]][g_seq[i]][g_seq[j + 1]];
}

double Ed5(int i, int j)
{
  return g_dangle5[g_seq[j]][g_seq[i]][g_seq[i - 1]];
}


double Ebi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > g_maxLoop)
    return 0.0;

  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] * g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]] * g_scale * g_scale;
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] * auPenalty(i, j) * auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] * pow(g_misc[12], log((double) loopSize2 / 30)) / g_scalen[loopSize2 - 30] * auPenalty(i, j) * auPenalty(ii, jj);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] * g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]] * g_scale * g_scale;
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] * auPenalty(i, j) * auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] * pow(g_misc[12], log((double) loopSize1 / 30)) / g_scalen[loopSize1 - 30] * auPenalty(i, j) * auPenalty(ii, jj);
    } /* Might leave these out (they are in the 2004 data tables, and we want 1999 for simplicity
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(g_seq[jj], g_seq[ii])][basePairIndex(g_seq[j], g_seq[i])][g_seq[jj + 1]][g_seq[ii - 1]][g_seq[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[i + 2]][g_seq[j - 2]];
  else if ((loopSize1 == 2 && loopSize2 == 3) ||
	   (loopSize1 == 3 && loopSize2 == 2))
    return g_tstacki23[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] *
      g_tstacki23[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]] / g_scalen[7];
      */
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] * pow(g_misc[12], log((double) (loopSize1 + loopSize2) / 30)) / g_scalen[loopSize1 + loopSize2 - 30];
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy *= g_tstacki[g_seq[i]][g_seq[j]][0][0];
	  loopEnergy *= g_tstacki[g_seq[jj]][g_seq[ii]][0][0];
	}
      else
	{
	  loopEnergy *= g_tstacki[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
	  loopEnergy *= g_tstacki[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
	}
      asPenalty = pow(g_misc[min3(4, loopSize1, loopSize2) - 1], abs(loopSize1 - loopSize2));
      if (asPenalty < g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy *= asPenalty;

      return loopEnergy;
    }
}
