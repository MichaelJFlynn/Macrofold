@echo off

gcc -Wall -std=c1x -g -o MacroFold StochasticSamples.h MacrofoldConsole.h PartitionFunction.h  AllowedPairs.h PairIterator.h EnergyModel.h EnergyFunctions.h RNA.h DataFile.h  EnergyModel.c RNA.c DataFile.c EnergyFunctions.c AllowedPairs.c PairIterator.c PartitionFunction.c MacrofoldConsole.c StochasticSamples.c main.c

gcc -Wall -std=c1x -g -o time_test StochasticSamples.h MacrofoldConsole.h PartitionFunction.h  AllowedPairs.h PairIterator.h EnergyModel.h EnergyFunctions.h RNA.h DataFile.h  EnergyModel.c RNA.c DataFile.c EnergyFunctions.c AllowedPairs.c PairIterator.c PartitionFunction.c MacrofoldConsole.c StochasticSamples.c time_test.c



