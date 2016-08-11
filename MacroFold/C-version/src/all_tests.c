// MacroFold includes
#include "EnergyFunctions.h"
#include "EnergyModel.h"
#include "RNA.h"

// Other includes
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>

// Some test sequences
char* seq1 = "AU";
char* seq2 = "AACCGAUCAUCAUCCCUCUCGAUGAUGAUGAUCCAUCAAAACACACCCUGUGUGUCGACUGCAGGGG";

// Other test variables

// Some useful functions
char intToBase(int num) {
  if (num == 0) 
    return 'A';
  else if (num == 1)
    return 'C';
  else if (num == 2)
    return 'G';
  else if (num == 3)
    return 'U';
  else
    return 'X';
}

 /* 
char* intToPair(int num) {
  if (num == 0) 
    return "AU";
  else if (num == 1)
    return "CG";
  else if (num == 2)
    return "GC";
  else if (num == 3)
    return "UA";
  else if (num == 4)
    return "GU";
  else if (num == 5)
    return "UG";
  else
    return "XX";
}
  
char* intToTwoBase(int num) {
  char first  = intToBase(num / 4);
  char second = intToBase(num % 4);
  char* third = (char *) malloc(3);
  strcpy(third, &first);
  strcat(third, &second);
  return third;
}
*/

void auPenaltyTest(char* strand) {
  RNA* my_strand = (RNA*) allocateRNA(strand);
  initializeEnergyModel(my_strand);
  // Print out the whole AU table for all possible combinations
  printf("AU Penalty Table:\n");
  for (int i=0; i < 16; ++i) {
    printf("AU Penalty for (%c, %c): %f\n", intToBase(i / 4), intToBase(i % 4), my_strand->energyModel->auPenalty[i/4][i%4]); 
  }
}

void hairpinTermTest(char* strand, int length) {
  // Print out the whole hairpinTest table for all possible combos
  printf("Hairpin Table\n");
  RNA* my_strand = (RNA*) allocateRNA(strand);
  initializeEnergyModel(my_strand);
  for (int i = 0; i < length; ++i) {
    for (int j = i+1; j < length; ++j) {
      if (hairpinTerm(my_strand, i, j) > 1e-6)
        printf("Energy for (%d, %d): %f\n", i, j, hairpinTerm(my_strand, i, j)); 
    }
  } 
}

void stackingTest() {
  char* string = malloc(5);
  string[4] = '\0';
  printf("Stacking Energies:\n");
  for (int i=0; i < 16; ++i) {
    for (int j=0; j < 16; ++j) {
      string[0] = intToBase(i/4);
      string[1] = intToBase(i%4);
      string[2] = intToBase(j%4);
      string[3] = intToBase(j/4);
      RNA* my_strand = (RNA*) allocateRNA(string);     
      initializeEnergyModel(my_strand);
      printf("Es for (%c, %c) then (%c, %c): %f\n", intToBase(i/4), intToBase(i%4), intToBase(j/4), intToBase(j%4), stackTerm(my_strand, 0, 3));
    }
  }
  free(string);
}

int main(int argc, char* argv[]) {

  //auPenaltyTest();
  hairpinTermTest(seq2, strlen(seq2));
  stackingTest();
}
