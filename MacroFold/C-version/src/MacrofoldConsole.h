
struct RNA;
typedef struct RNA RNA;

typedef struct MacrofoldConsole {
  RNA** strands;
  int numStrands;
} MacrofoldConsole;

void startConsole(MacrofoldConsole* mc);
int execLine(MacrofoldConsole* mc);
MacrofoldConsole* allocateMacrofoldConsole();

void freeMacrofoldConsole(MacrofoldConsole* mc);
