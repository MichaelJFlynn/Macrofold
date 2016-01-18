@echo off

gcc -Wall -o MacroFold EnergyModel.h EnergyFunctions.h RNA.h DataFile.h EnergyModel.c RNA.c DataFile.c EnergyFunctions.c main.c