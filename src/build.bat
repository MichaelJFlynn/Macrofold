@echo off

gcc -Wall -o MacroFold PartitionFunction.h  AllowedPairs.h PairIterator.h EnergyModel.h EnergyFunctions.h RNA.h DataFile.h EnergyModel.c RNA.c DataFile.c EnergyFunctions.c AllowedPairs.c PairIterator.c PartitionFunction.c main.c