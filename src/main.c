// Written by Mike Flynn

#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include "EnergyFunctions.h"
#include "PartitionFunction.h"
#include "AllowedPairs.h"
#include "MacrofoldConsole.h"
#include "StochasticSamples.h"
#include <stdio.h>
#include <math.h>

/* 

Notes:
- new scaling factor (for now .34 kcal/mol*base / kT)
- markham: 2d array in 1d array?
- Zuker: markham parraellized


Official list of tests:
- compare free energy to Unafold
- compare P[i][j] to unafold
- simple structures (small hairpins) does it output hairpin
- compare longer sequences for multibranch rules
- compare to cannonical RNA sequences:
---- TRNAs, 5s Ribosomal RNA -> RNAStructure examples dir, 09wkj Nestor sequences dir


1: implement AU
2: implement tstackh
3: implement gscale
4: implement multiloop energies (can find in nndb)
5: unthorough verify
5.5?: if 1x2 mismatches differ: implement 1x2 mismatches and 2x2 mismatches
5: timing tests, random sequence stuff, go by 100s or 1000, out to 4000
6: dynamic thresholding
7: test again?
8: Get paper together - what do we need?
9: thorough verify against unafold, test above
 */


int main(int argc, char* argv[]) {
  // tests 
  RNA* strand = readSequenceFile("sample.seq");

  fillZbZ1Z2(strand);
  fillZ(strand);
  fillExtendedZbZ1Z2(strand);
  fillP(strand);

  //printf("%g\n", getFreeEnergy(strand));

  /* MacrofoldConsole* mc =  allocateMacrofoldConsole(); */
  /* startConsole(mc); */


  sample(strand, 1000);
  int i, j;
  for(i = 0; i < 1000; i++) {
    for(j = 0; j < strand->length; j++) {
      printf("%d ", strand->samples->structures[i][j]);
    }
    printf("\n");    
  }



  /* free(mc); */
  /* scanf("%s", scan_string); */
  /* printf("You entered: %s adsf\n", scan_string); */

  /* int i,j; */
  /* double sum = 0; */
  /* for(i = 0; i < strand->length; i++) { */
  /*   printf("%11d\t", i); */
  /* } */
  /* printf("\n"); // */
  /* for(i = 0; i < strand->length; i++) { */
  /*   for(j = 0; j < strand->length; j++) { */
  /*     if(i < j) { */
  /* 	printf("%.4e\t", strand->partitionFunction->P[i][j]); */
  /*     } else { */
  /* 	printf("           \t"); */
  /*     } */
  /*   } */
  /*   printf("\n"); */
  /* } */

  //i = 0;
  //j = 6; 
  //double** Zb = strand->partitionFunction->Zb; 
  //printf("Zb[%d, %d] = %g\n", i, j, hairpinTerm(strand, i, j));


  /* for(i = 0; i < strand->length; i++) { */
  /*   for(j = 0; j < strand->length; j++) { */
  /*     //      double prob = strand->partitionFunction->P[i][j]; */
  /*     if(strand->partitionFunction->Zb[i][j] > 0) */
  /* 	printf("P[%d][%d] = (%g)(%g)/(%g) =  %g\n", i, j, */
  /* 	       strand->partitionFunction->Zb[i][j], */
  /* 	       strand->partitionFunction->Zb[j][i + strand->length], */
  /* 	       strand->partitionFunction->Z[0][strand->length - 1], */
  /* 	       strand->partitionFunction->P[i][j]); */
  /*   } */
  /* } */

  

  /* double** Zb = strand->partitionFunction->Zb; */
  /* double** Z2 = strand->partitionFunction->Z2; */
  /* double** Z = strand->partitionFunction->Z; */
  /* double* scale = strand->energyModel->scale; */
  /* int len = strand->length; */
  /* double multiA = 1, multiB = 1, multiC = 1, au = 1; */
  /* i = 26; */
  /* j = 10 + len; */
  /* printf("Zb[%d][%d] %g = %g + %g\n", i, j, Zb[i][j], stackTerm(strand, i, j - len) * Zb[i + 1][j - 1],  bulgeInternalTerm(strand, i, j - len)); */
  /* printf("Zb[%d][%d] +=  %g\n", i, j, multiA * multiC * au * (Z2[i+1][j-1] */
  /* 							      + ed5(strand, j - len, i) * multiB *Z2[i+1][j-2] */

  /* 							      + etstackm(strand, j - len,i) * multiB * multiB * Z2[i+2][j-2])); */
  /* printf("Zb[%d][%d] += %g * (%g + %g) * (%g + %g)/ %g\n", i, j,  */
  /* 	 au, Z[i+1][len - 1], 1/scale[len - 1 - i], */
  /* 	 Z[0][j - 1 - len], 1/scale[j-len], scale[2]); */

  /* printf("Zb[%d][%d] += %g * %g * (%g + %g) * (%g + %g)/%g\n", i, j, au, ed5(strand, j - len, i),Z[i + 1][len - 1], 1/scale[len - 1 - i], */
  /* 	 Z[0][j - 2 - len], 1/scale[j-len-1], scale[2]);  */

  /* printf("Zb[%d][%d] += %g * %g * (%g + %g) * (%g + %g)/%g\n", i, j, */
  /* 	 au, ed3(strand, j - len, i), Zb[i + 2][len - 1], 1/scale[len - 2 - i], */
  /* 	 Zb[0][j - 1 - len], 1/scale[j - len], scale[2]); */

  /* printf("Zb[%d][%d] += %g * %g * (%g + %g) * (%g + %g)/%g\n", i,j, */
  /* 	 au, etstackm(strand, j - len, i), Zb[i+2][len-1], 1/scale[len - 2 - i],  */
  /* 	 Zb[0][j - 2 - len], 1/scale[j-len-1],scale[2]); */
  /* printf("P(i,j) = %g\n", Zb[10][26] * Zb[26][50] / Z[0][strand->length - 1]); */
  /* printf("%g\n",  Zb[9][27]); */
  /* printf("%g\n",  stackTerm(strand, 26, 10)); */
  /* int* nSeq = strand->intSequence; */
  /* printf("%g\n", strand->energyModel->stack[nSeq[26]][nSeq[10]][nSeq[26 + 1]][nSeq[10 - 1]]); */
  /* printf("%g\n", strand->energyModel->stack[nSeq[10]][nSeq[26]][nSeq[10 + 1]][nSeq[26 - 1]]); */
  
  /*
    Reporting Free Energies
    Finding Pijs
    Stochastic traceback
   */

  freeRNA(strand);
  return 0;
}
