# Ternary LWE Parameter #

## Generate Samples ##
```
% time ../../bin/lwesampler -n80 -m120 -s1.0 -q4093 -3
../../bin/lwesampler -n80 -m120 -s1.0 -q4093 -3  0.02s user 0.00s system 98% cpu 0.018 total
```

## Reduction ##
```
% time ../../bin/reduction -b2
---- LLL_RR status ----
elapsed time: 10:50, stage: 121, rank: 120, swaps: 332116
log of prod of lengths: 917.327
dumping to BKZ_dump.dat...
---- BKZ_RR status ----
elapsed time: 10:50, enum time: 0, iter: 119
triv: 0, nontriv: 0, no ops: 119, rank: 120, swaps: 332116
log of prod of lengths: 917.327
dumping to BKZ_dump.dat...
../../bin/reduction -b2  649.43s user 0.24s system 99% cpu 10:49.89 total

% cp samples_vector.dat reduced_vector.dat
```
Note that you have to copy the samples_vector file to reduced_vector file, because
the enumeration program looks for the matrix/vector in files with the same prefix.

## Enumeration ##
```
% time ../../bin/enumeration -n80 -s1.0 -q4093 -b2 -eln -i"reduced" -p --factor=1.0
running parallel implementation
breadth-first traversal till k = 111
Liu Nguyen			[done]
../../bin/enumeration -n80 -s1.0 -q4093 -b2 -eln -i"reduced" -p --factor=1.0  3.07s user 0.02s system 103% cpu 3.000 total
```

## Test Solution ##
```
% diff solution_vector.dat samples_solution.dat
```
For the correct solution, this results in no difference.
