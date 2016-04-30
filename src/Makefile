all:	
	gcc -std=c11 -o MacroFold EnergyModel.c RNA.c DataFile.c EnergyFunctions.c AllowedPairs.c PairIterator.c PartitionFunction.c MacrofoldConsole.c StochasticSamples.c main.c

test_pij:
	gcc -std=c11 -g -o test_pij EnergyModel.c RNA.c DataFile.c EnergyFunctions.c AllowedPairs.c PairIterator.c PartitionFunction.c MacrofoldConsole.c StochasticSamples.c test_pij.c
