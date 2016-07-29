// MacroFold includes
#include "EnergyFunctions.h"
#include "EnergyModel.h"
#include "RNA.h"

// UNAFold includes
//#include "../../unafold-3.9/src/energy.h"

// Other includes
#include <stdio.h>
#include <math.h>
#include <time.h>

// Some test sequences
char* seq1 = "AA";
/*char* seq2 = "AC";
char* seq3 = "AG";
char* seq4 = "AU";
char* seq5 = "CA";
char* seq6 = "CC";
char* seq7 = "CG";
char* seq8 = "CU";
char* seq9 = "GA";
char* seq10 = "GC";
char* seq11 = "GG";
char* seq12 = "GU";
char* seq13 = "UA";
char* seq14 = "UC";
char* seq15 = "UG";
char* seq16 = "UU";
*/
void auPenaltyTest(char* seq) {
  // We start with MacroFold

  // Load in a strand
  RNA* my_strand = (RNA*) allocateRNA(seq);

  // Initialize the Energy Model for that strand
  initializeEnergyModel(my_strand);

  // Save the auPenalty for that strand
  printf("%f",my_strand->energyModel->auPenalty[seq[0]][seq[1]]);

  // We do the same with UNAFold
  //double aup_unafold = auPenalty(seq[0], seq[1]);

  // Check if energies are the same
  /*if (aup_macrofold == aup_unafold)
    exit(0);
  else {
    printf("AU Penalty Test Failed <Details>\n");
    exit(-1);
    */
 
}

int main(int argc, char* argv[]) {

  auPenaltyTest(seq1);

}
