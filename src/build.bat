@echo off

gcc -Wall -o MacroFold MacrofoldConsole.h PartitionFunction.h  AllowedPairs.h PairIterator.h EnergyModel.h EnergyFunctions.h RNA.h DataFile.h  EnergyModel.c RNA.c DataFile.c EnergyFunctions.c AllowedPairs.c PairIterator.c PartitionFunction.c MacrofoldConsole.c main.c


