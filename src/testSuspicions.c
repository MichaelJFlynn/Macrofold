#include <stdio.h>

typedef struct { 
  char* str;
} holder;

typedef struct {
  char stack[4][4][4][4];
} arrayHolder;

holder funcOnHolder(holder hold) {
  char* word = "You are transformed!";
  hold.str = word;
  return hold;
}

int main(int argc, char** argv) {
  holder hold;
  hold.str = "I am not transformed!";
  hold = funcOnHolder(hold);
  arrayHolder ah;
  ah.stack[1][2][3][0] = 'c';
  printf("%s\n%lu\n", hold.str, sizeof(arrayHolder));
}
