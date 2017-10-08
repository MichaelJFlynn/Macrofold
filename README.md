MacroFold
Michael Flynn & Daniel Aalberts

This implementation is meant to demonstrate algorithmic
improvements. The 1999 free energy rules are implemented.

To Build
-------------

On Unix based systems (with gcc installed):

```
cd src 
make all
```

On Windows systems with gcc installed:

``` 
cd src 
./build.bat
```


Example Usage
----------------

From `src/time_test.c`:

```c
#include "RNA.h"

int main(int argc, char* argv[]) {

  RNA* strand = readSequenceFile(argv[1]);
  computePartitionFunction(strand);
}
```