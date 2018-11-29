- This is an initial submission 
- One note: 

"1 NOTE
checking compiled code ... NOTE
File ‘FastLZeroSpikeInference/libs/FastLZeroSpikeInference.so’:
  Found ‘_printf’, possibly from ‘printf’ (C)
    Objects: ‘ARFPOP.o’, ‘funPieceListLog.o’
  Found ‘_puts’, possibly from ‘printf’ (C), ‘puts’ (C)
    Objects: ‘ARFPOP.o’, ‘funPieceListLog.o’

Compiled code should not call entry points which might terminate R nor
write to stdout/stderr instead of to the console, nor the system RNG."

This note cannot easily be resolved, because the compiled library FastLZeroSpikeInference.so must be compilable on systems *without* R headers (which would be required to replace all printf commands with Rprintf commands). This is because this R package is simply a wrapper to the compiled library FastLZeroSpikeInference.so. Other languages, i.e. python, need to be able to call this library and so introducing additional libraries is highly undiserable. 

