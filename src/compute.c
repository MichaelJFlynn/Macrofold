// Written by Mike Flynn, 2015

/*
TODO: 
*fields that strand needs:
** length
** sequence (obviously)
 */

#include "EnergyFunctions.h"

#define TURN 3

#ifdef DEBUG


#else

#endif


RNA computePartitionFunction(RNA strand) {
  int i, j, k, index;
  
  double multiA = strand.energyModel.multiA;
  double multiB = strand.energyModel.multiB;
  double multiC = strand.energyModel.multiC;
  double scale = strand.energyModel.scalingTerm;

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

  for (j = 2; j <= strand.length; ++j)
    for (i = j - TURN - 1; i >= 1; --i)
      {
	double au = auPenalty(i, j);
	double Qmulti = 0;
	double Zmulti = 0; 

	if(Zb(i, j) != 0.0) {
	  Zb(i, j) = HairpinTerm(i, j) + StackTerm(i,j) * Zb(i + 1, j - 1) + BulgeInternalTerm(i, j);	
	  Zmulti += multiA * multiC * au * (Z2(i+1, j-1) 
					     + Ed5(j,i) * multiB *Z2(i+1, j-2)
					     + Ed3(j,i) * multiB * Z2(i+2, j-1) 
					     + Etstackm(j,i) * multiB * multiB * Z2(i+2, j-2));
	  Zb(i, j) += Zmulti;
	}

	Z2(i, j) = multiB * Z2(i+1, j) / scale;
	Z1(i, j) = multiB * Z1(i+1, j) / scale;
	for(index=0; pfc.ij[i][index] != 0 && pfc.ij[i][index] <= j; ++index) {
	  k = pfc.ij[i][index];
	  Z1(i,j) += multiC * auPenalty(i, k) * Zb(i, k) * g_bscalen[j-k];	  
	  if(k < j) { 
	    Z1(i,j) += multiC * auPenalty(i, k) * Zb(i, k) * Z1(k+1, j); 
	    Z1(i,j) += multiC * multiB * auPenalty(i, k) * Zb(i, k) * g_bscalen[j-(k+1)] * Ed3(i, k);
	    Z2(i,j) += multiC * auPenalty(i, k) * Zb(i, k) * Z1(k+1, j);	  
	    if(k < j - 1) {
	      Z1(i,j) += multiC * multiB * auPenalty(i, k) * Zb(i, k) * Z1(k+2, j) * Ed3(i, k);
	      Z2(i,j) += multiC * multiB * auPenalty(i, k) * Zb(i, k) * Z1(k+2, j) * Ed3(i, k);
	    }
	  }
	}

	for(index=0; pfc.ij[i+1][index] != 0 && pfc.ij[i+1][index] <= j; ++index) {
	  k = pfc.ij[i+1][index];
	  Z1(i,j) += multiC * multiB * auPenalty(i+1, k) * Zb(i+1, k) * g_bscalen[j-k] * Ed5(i+1, k);
	  if(k < j) {
	    Z1(i,j) += multiC * multiB * auPenalty(i+1, k) * Zb(i+1, k) * Z1(k+1, j) * Ed5(i+1, k);
	    Z1(i,j) += multiC * multiB * multiB * auPenalty(i+1, k) * Zb(i+1, k) * g_bscalen[j- (k+1)] * Etstackm(i+1, k);
	    Z2(i,j) += multiC * multiB * auPenalty(i+1, k) * Zb(i+1, k) * Z1(k+1, j) * Ed5(i+1, k);
	    if(k < j - 1) {
	      Z1(i,j) += multiC * multiB * multiB * auPenalty(i+1, k) * Zb(i+1, k) * Z1(k+2, j) * Etstackm(i+1, k);	  
	      Z2(i,j) += multiC * multiB * multiB * auPenalty(i+1, k) * Zb(i+1, k) * Z1(k+2, j) * Etstackm(i+1, k);
	    }
	  }
	}

      }
}


void computeQ53FAST(double* q5, double* q3) {
  int i, j, index, k;
  Z(1,1) = 0.0;
  Z(strand.length, strand.length) = 0.0;

  for (j = 2; j <= strand.length; ++j) { 

    Z(1,j) = Z(1, j - 1) / scale;
    for(index=0; pfc.ji[j][index] != 0 && pfc.ji[j][index] >= 1; ++index) {
      k = pfc.ji[j][index];
      Z(1,j) += auPenalty(k, j) * Zb(k, j) / scalen[k - 1];
      if(k > 1) {
	Z(1,j) += auPenalty(k, j) * Z(1, k-1) * Zb(k, j); 
	Z(1,j) += auPenalty(k, j) * Zb(k, j) / scalen[k-2] * Ed5(k, j);
	if(k > 2) {
	  Z(1,j) += auPenalty(k, j) * Z(1, k-2) * Zb(k, j) * Ed5(k, j);
	}
      }
    }
      
    for(index=0; pfc.ji[j-1][index] != 0 && pfc.ji[j-1][index] >= 1; ++index) {
      k = pfc.ji[j-1][index];
      Z(1,j) += auPenalty(k, j-1) * Zb(k, j-1) * Ed3(k, j-1) / scalen[k-1];
      if(k > 1) {
	Z(1,j) += auPenalty(k, j-1) * Z(1, k-1) * Zb(k, j-1) * Ed3(k, j-1);
	Z(1,j) += auPenalty(k, j-1) * Zb(k, j-1) * Etstackm(k, j-1) / scalen[k-2];	  
	if(k > 2) {
	  Z(1,j) += auPenalty(k, j-1) * Z(1, k-2) * Zb(k, j-1) * Etstackm(k, j-1);	  
	}
      }
    
    }
  }

  for (i = strand.length - 1; i >= 1; --i) {
    Z(i, strand.length) = Z(i + 1, strand.length) / scale;

    for(index=0; pfc.ij[i][index] != 0 && pfc.ij[i][index] <= strand.length; ++index) {
      k = pfc.ij[i][index];
      Z(i,strand.length) += auPenalty(i, k) * Zb(i, k)/scalen[strand.length - k];
      if(k < strand.length) {
	Z(i,strand.length) += auPenalty(i, k) * Zb(i, k) * Z(k+1, strand.length); 
	Z(i, strand.length) += auPenalty(i, k) * Zb(i, k) / scalen[strand.length - k - 1] * Ed3(i,k);
	if(k < strand.length-1) {
	  Z(i,strand.length) += auPenalty(i, k) * Zb(i, k) * Z(k+2, strand.length) * Ed3(i, k);
	}
      }
    }

    for(index=0; pfc.ij[i+1][index] != 0 && pfc.ij[i+1][index] <= strand.length; ++index) {
      k = pfc.ij[i+1][index];
      Z(i, strand.length) += auPenalty(i+1, k) * Zb(i+1, k) / scalen[strand.length - k] * Ed5(i+1, k);
      if(k < strand.length) {
	Z(i,strand.length) += auPenalty(i+1, k) * Zb(i+1, k) * Z(k+1, strand.length) * Ed5(i+1, k); 
	Z(i,strand.length) += auPenalty(i+1, k) * Zb(i+1, k) / scalen[strand.length - k - 1] * Etstackm(i + 1, k);
	if(k < strand.length-1) {
	  Z(i,strand.length) += auPenalty(i+1, k) * Zb(i+1, k) * Z(k+2, strand.length) * Etstackm(i + 1, k);
	}
      }
    }
  }
}


double* Z(int i, int j) {

}

double* Z1(int i, int j) {

}

double* Zb(int i, int j) {

}

double* Z2(int i, int j) {

}
