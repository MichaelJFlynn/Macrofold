#ifndef ENERGYFUNCTIONS_H
#define ENERGYFUNCTIONS_H

#define TURN 3

// forward declaration of RNA, for the compiler
struct RNA;
typedef struct RNA RNA;

double hairpinTerm(RNA* strand, int i, int j);
double bulgeInternalTerm(RNA* strand, int i, int j);
double bulgeInternalTermPeriodic(RNA* strand, int i, int j);
double stackTerm(RNA* strand, int i, int j);
double etstackm(RNA* strand, int i, int j);
double ed3(RNA* strand, int i, int j);
double ed5(RNA* strand, int i, int j);
double ebi(RNA* strand, int i, int j, int ii, int jj);


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
 * Implement scaling parameter (unafold has estimateScalingParameter())
 */

#endif
