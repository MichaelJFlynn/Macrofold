// Written by Mike Flynn

#include "EnergyModel.h"
#include "RNA.h"
#include "DataFile.h"
#include "EnergyFunctions.h"
#include "PartitionFunction.h"
#include "AllowedPairs.h"
#include <stdio.h>
#include <math.h>


int main(int argc, char* argv[]) {
  // tests 
  RNA* strand = readSequenceFile("sample.seq");

  fillZbZ1Z2(strand);
  fillZ(strand);
  fillExtendedZbZ1Z2(strand);
  fillP(strand);

  printf("%g\n", getFreeEnergy(strand));

  int i,j;
  //double sum = 0;
  /* for(i = 0; i < strand->length; i++) { */
  /*   printf("%11d\t", i); */
  /* } */
  /* printf("\n"); // */
  /* for(i = 0; i < strand->length; i++) { */
  /*   for(j = 0; j < strand->length; j++) { */
  /*     if(i < j) { */
  /* 	printf("%.4e\t", strand->partitionFunction->Zb[i][j]); */
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


  for(i = 0; i < strand->length; i++) {
    for(j = 0; j < strand->length; j++) {
      double prob = strand->partitionFunction->P[i][j];
      if(prob > .0001)
  	printf("P[%d][%d] = %g\n", i, j, prob);
    }
  }
  /* double** Zb = strand->partitionFunction->Zb; */
  /* double** Z2 = strand->partitionFunction->Z2; */
  /* double multiA = 1, multiB = 1, multiC = 1, au = 1; */
  /* i = 0; */
  /* j = 6;  */
  /* printf("Zb[%d][%d] = %g\n", i, j, hairpinTerm(strand, i, j) + stackTerm(strand, i,j) * Zb[i + 1][j - 1] + bulgeInternalTerm(strand, i, j));	 */
  /* printf("Zb[%d][%d] +=  %g\n", i, j, multiA * multiC * au * (Z2[i+1][j-1]  */
  /* 				      + ed5(strand, j, i) * multiB *Z2[i+1][j-2] */
  /* 				      + ed3(strand, j,i) * multiB * Z2[i+2][j-1]  */
  /* 							    + etstackm(strand, j,i) * multiB * multiB * Z2[i+2][j-2]));	 */
  

  /*
    Reporting Free Energies
    Finding Pijs
    Stochastic traceback
   */

  freeRNA(strand);
  return 0;
}
