#ifndef DATAFILE_H
#define DATAFILE_H

typedef struct {
  char** tokens;
  int length;
} Line;

typedef struct {
  Line* header;
  Line** data;
  int nrow;
} DataFile;


DataFile* readCSV(char* filepath);

char* get(DataFile* data, int i, int j);
Line* parseLine(char* string, char* tokens);
/* inline int nrow(DataFile data) { */
/*   return data.nrow; */
/* };  */

/* inline int nrowp(DataFile* data) { */
/*   return data->nrow; */
/* }; */

/* inline int ncol(DataFile data) { */
/*   return data.header->length; */
/* }; */

/* inline int ncolp(DataFile* data) { */
/*   return data->header->length; */
/* }; */

DataFile* constructDataFile(Line** lineHolder, int nrow);

void freeDataFile(DataFile* data);
void freeLine(Line* line);
#endif
