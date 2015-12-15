


double dangleEnergies3[4][4][4];
double dangleEnthalpies3[5][5][6];
double dangleEnergies5[4][4][4];
double dangleEnthalpies5[5][5][6];
double stackEnergies[4][4][4][4];
double stackEnthalpies[5][5][5][5];
double interiorLoopEnergies[30];
double bulgeLoopEnergies[30];
double hairpinLoopEnergies[30];
double interiorLoopEnthalpies[30];
double bulgeLoopEnthalpies[30];
double hairpinLoopEnthalpies[30];
double sint2Energies[6][6][4][4];
double sint2Enthalpies[7][7][5][5];
double asint1x2Energies[6][6][4][4][4];
double asint1x2Enthalpies[7][7][5][5][5];
double sint4Energies[6][6][4][4][4][4];
double sint4Enthalpies[7][7][5][5][5][5];
double tstackhEnergies[4][4][4][4];
double tstackhEnthalpies[5][5][5][5];
double tstackiEnergies[4][4][4][4];
double tstackiEnthalpies[5][5][5][5];
double tstacki23Energies[4][4][4][4];
double tstacki23Enthalpies[5][5][5][5];
double tstackmEnergies[4][4][4][4];
double tstackmEnthalpies[5][5][6][6];
double tstackeEnergies[4][4][4][4];
double tstackeEnthalpies[5][5][6][6];
struct triloopE* triloopEnergies;
struct triloopE* triloopEnthalpies;
struct tloopE* tloopEnergies;
struct tloopE* tloopEnthalpies;
struct hexaloopE* hexaloopEnergies;
struct hexaloopE* hexaloopEnthalpies;
double multiEnergies[3];
double multiEnthalpies[3];
double miscEnergies[13];
double miscEnthalpies[13];


void loadData() {
  
  loadStack(stackEnergies, stackEnthalpies, NA, saltCorrection);
  symmetryCheckStack(stackEnergies, "energy");
  /* symmetryCheckStack(stackEnthalpies, "enthalpy"); */
  loadDangle(dangleEnergies3, dangleEnthalpies3, dangleEnergies5, dangleEnthalpies5, NA, saltCorrection);
  loadLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, NA, saltCorrection);
  loadSint2(sint2Energies, sint2Enthalpies, NA, saltCorrection);
  symmetryCheckSint2(sint2Energies, "energy");
  /* symmetryCheckSint2(sint2Enthalpies, "enthalpy"); */
  loadAsint1x2(asint1x2Energies, asint1x2Enthalpies, NA, saltCorrection);
  loadSint4(sint4Energies, sint4Enthalpies, NA, saltCorrection);
  symmetryCheckSint4(sint4Energies, "energy");
  /* symmetryCheckSint4(sint4Enthalpies, "enthalpy"); */
  loadTstackh(tstackhEnergies, tstackhEnthalpies, NA);
  loadTstacki(tstackiEnergies, tstackiEnthalpies, NA);
  loadTstacki23(tstacki23Energies, tstacki23Enthalpies, NA);
  loadTstackm(tstackmEnergies, tstackmEnthalpies, NA, saltCorrection);
  loadTstacke(tstackeEnergies, tstackeEnthalpies, NA, saltCorrection);
  loadTriloop(&triloopEnergies, &triloopEnthalpies, &numTriloops, NA);
  g_triloop = (struct triloop*) xcalloc(numTriloops, sizeof(struct triloop));
  loadTloop(&tloopEnergies, &tloopEnthalpies, &numTloops, NA);
  g_tloop = (struct tloop*) xcalloc(numTloops, sizeof(struct tloop));
  loadHexaloop(&hexaloopEnergies, &hexaloopEnthalpies, &numHexaloops, NA);
  g_hexaloop = (struct hexaloop*) xcalloc(numHexaloops, sizeof(struct hexaloop));
  loadMulti(multiEnergies, multiEnthalpies, NA);
  loadMisc(miscEnergies, miscEnthalpies, NA);
  
  combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
  combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
  combineLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, tRatio, g_hairpinLoop, g_interiorLoop, g_bulgeLoop);
  combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
  combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
  combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
  combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
  combineTstack(tstacki23Energies, tstacki23Enthalpies, tRatio, g_tstacki23);
  combineTstack(tstackhEnergies, tstackhEnthalpies, tRatio, g_tstackh);
  
  combineTstack2(tstackmEnergies, tstackmEnthalpies, tRatio, g_tstackm);
  combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
  
  combineMulti(multiEnergies, multiEnthalpies, tRatio, g_multi);
  combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
  combineTriloop(triloopEnergies, triloopEnthalpies, tRatio, g_triloop, numTriloops);
  combineTloop(tloopEnergies, tloopEnthalpies, tRatio, g_tloop, numTloops);
  combineHexaloop(hexaloopEnergies, hexaloopEnthalpies, tRatio, g_hexaloop, numHexaloops);
  
  if (g_scale < 1.0)
    {
      g_scale = exp(-estimateScale(g_stack) / RT);
      if (g_scale < 1.0)
	g_scale = 1.0;
    }
  free(g_scalen);
  g_scalen = xcalloc(g_len + 1, sizeof(double));
  g_scalen[0] = 1.0;
  for (i = 1; i <= g_len; ++i)
    g_scalen[i] = g_scalen[i - 1] * g_scale;
  
  calculateStack(g_stack, tRatio, g_scale);

  calculateZipDangle(g_dangle3, g_dangle5, tRatio, g_scale);
  calculateDangle(g_dangle3, g_dangle5, tRatio, g_scale);
  calculateLoop(g_hairpinLoop, g_interiorLoop, g_bulgeLoop, tRatio, g_scale);
  calculateSint2(g_sint2, tRatio, g_scale);
  calculateAsint1x2(g_asint1x2, tRatio, g_scale);
  calculateSint4(g_sint4, tRatio, g_scale);
  calculateStack(g_tstacki, tRatio, 1.0);
  calculateStack(g_tstacki23, tRatio, 1.0);
  calculateStack(g_tstackh, tRatio, 1.0);
  calculateStack2(g_tstackm, tRatio, g_scale);
  calculateStack2(g_tstacke, tRatio, g_scale);
  calculateMulti(g_multi, tRatio, g_scale);
  calculateMisc(g_misc, tRatio);
  calculateTriloop(g_triloop, numTriloops, tRatio);
  calculateTloop(g_tloop, numTloops, tRatio);
  calculateHexaloop(g_hexaloop, numHexaloops, tRatio);
}
