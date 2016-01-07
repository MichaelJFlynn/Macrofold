#include "DataFile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Util files for parsing data
   
   Not meant for general usage, rather the specific parsing of our
   datafiles. Assumptions:
   ---- data files have a header in the first line, seperated by 
   ----- commas
   ---- all data file lines are < 1000 characters long
 */

#define BUFSIZE 1000

Line* parseLine(char* line, char* tokens) {
  char** parsedLine;
  char* tokenHolder[BUFSIZE];
  char* token;
  int i, numTokens = 0;
  Line* result = (Line*) malloc(sizeof(Line));

  // calculate the number of tokens with the disposable data
  token = strtok(line, tokens);
  tokenHolder[0] = (char*) malloc(strlen(token) + 1);
  memcpy(tokenHolder[0], token, strlen(token) + 1);
  numTokens++;
  while((token = strtok(NULL, tokens))) {
    tokenHolder[numTokens] = (char*) malloc(strlen(token) + 1);
    memcpy(tokenHolder[numTokens], token, strlen(token) + 1);
    numTokens++;
  }

  // now actually create the char**
  parsedLine = (char**) malloc(numTokens * sizeof(char*));
  for(i = 0; i < numTokens; i++) {
    parsedLine[i] = tokenHolder[i];
  }

  result->tokens = parsedLine;
  result->length = numTokens;
  return result;
}

DataFile* readCSV(char* filepath) {
  FILE* fp;
  char buf[BUFSIZE];
  Line* lineHolder[BUFSIZE];
  int i = 0;

  // open file
  fp = fopen(filepath, "r");
  if(!fp) {
    printf("Could not open file: %s\n", filepath);
    exit(1);
  }

  while(fgets(buf, BUFSIZE, fp)) {    
    lineHolder[i] = parseLine(buf, ",\n");
    i++;
  }

  // close file, reading done
  fclose(fp);

  return constructDataFile(lineHolder, i);
}

char* get(DataFile* df, int i, int j) {
  if(i >= df->nrow) {
    printf("illegal access, nrows: %d, i = %d.\n", df->nrow, i);
    exit(1); 
  }
  
  if(j >= df->header->length) {
    printf("illegal access, ncols: %d, j = %d.\n", df->header->length, j);
    exit(1);
  }

  return df->data[i]->tokens[j];
}

// memory could be consildated here, at least defragmented, at this
// step, if memory ever becomes an issue. 
DataFile* constructDataFile(Line** lines, int nrow) {
  DataFile* df = (DataFile*) malloc(sizeof(DataFile));
  df->data = (Line**) malloc((nrow- 1) * sizeof(Line*));
  int i;

  df->header = lines[0];
  df->nrow = nrow - 1;
  for(i = 1; i < nrow; i++) {
    if(lines[i]->length != df->header->length) {
      printf("Inconsistant column number: %d, %d\n", lines[i]->length, df->header->length);
    }
    df->data[i-1] = lines[i];
  }
  return df;
}

void freeDataFile(DataFile* dataFile) {
  int i, j;

  for(i = 0; i < dataFile->nrow; i++) {    
    for(j = 0; j < dataFile->header->length; j++) {
      free(dataFile->data[i]->tokens[j]);
    }
    free(dataFile->data[i]->tokens);
    free(dataFile->data[i]);
  }

  for(i = 0; i < dataFile->header->length; i++) {
    free(dataFile->header->tokens[i]);
  }
  free(dataFile->header->tokens);
  free(dataFile->header);

  free(dataFile->data);
 
  free(dataFile);
}
