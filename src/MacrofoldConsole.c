#include "RNA.h"
#include "MacrofoldConsole.h"
#include "DataFile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int MAXSTRANDS = 1000;

void startConsole(MacrofoldConsole* mc) {
  
  printf(" _  _   __    ___  ____   __  ____  __   __    ____ \n"
	 "( \\/ ) / _\\  / __)(  _ \\ /  \\(  __)/  \\ (  )  (    \\\n"
	 "/ \\/ \\/    \\( (__  )   /(  O )) _)(  O )/ (_/\\ ) D (\n"
	 "\\_)(_/\\_/\\_/ \\___)(__\\_) \\__/(__)  \\__/ \\____/(____/\n" 
	 "\nWelcome to the Macrofold Console!\n"
	 "Written by Mike Flynn <mflynn210@gmail.com> 2015-2016\n"
	 "Advisor: Daniel P. Aalberts <aalberts@williams.edu>\n"
	 "\nType 'help' for a list of commands.\n\n"
	 );
  do {
    printf("> ");
  } while(execLine(mc));
  
}

void printHelp() {
  printf("%20s\t%-50s\n\n", "---- Command", "Description ----");
  printf("%20s\t%-50s\n", "help", "Prints help information");
  printf("%20s\t%-50s\n", "exit", "Exits console");
}


int execLine(MacrofoldConsole* mc) {
  char buff[500];
  int finishCode = 1;
  fflush(stdout);
  fgets(buff, 500, stdin);
  
  Line* line = parseLine(buff, " \n");
  
  if(line->length == 0) {
    return 1;
  }
  char* command = line->tokens[0];

  if(strcmp(command, "exit") == 0) {
    finishCode = 0;
  } else if(strcmp(command, "help") == 0) {
    printHelp();
  } else {
    printf("Unknown command: \"%s\"\n", command);
  }

  freeLine(line);
  return finishCode;
}


MacrofoldConsole* allocateMacrofoldConsole() {
  MacrofoldConsole* newMC = (MacrofoldConsole *) malloc(sizeof(MacrofoldConsole));
  newMC->strands = malloc(MAXSTRANDS * sizeof(RNA*));
  return newMC;
}

void freeMacrofoldConsole(MacrofoldConsole* mc) { 
  free(mc->strands); 
  free(mc);
}
