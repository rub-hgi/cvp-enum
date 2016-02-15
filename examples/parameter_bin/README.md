# Binary LWE Parameter #

## Generate Samples ##
```
% time ../../bin/lwesampler -n80 -m120 -s1.0 -q4093 -2 
../../bin/lwesampler -n80 -m120 -s1.0 -q4093 -2  0.02s user 0.00s system 96% cpu 0.019 total
```

## Reduction ##
```
% time ../../bin/reduction -b2
---- LLL_RR status ----
elapsed time: 11:02, stage: 121, rank: 120, swaps: 333382
log of prod of lengths: 916.975
dumping to BKZ_dump.dat...
---- BKZ_RR status ----
elapsed time: 11:02, enum time: 0, iter: 119
triv: 0, nontriv: 0, no ops: 119, rank: 120, swaps: 333382
log of prod of lengths: 916.975
dumping to BKZ_dump.dat...
../../bin/reduction -b2  661.47s user 0.24s system 99% cpu 11:01.93 total

% cp samples_vector.dat reduced_vector.dat
```
Note that you have to copy the samples_vector file to reduced_vector file, because
the enumeration program looks for the matrix/vector in files with the same prefix.

## Enumeration ##
```
% time ../../bin/enumeration -n80 -s1.0 -q4093 -b2 -eln -i"reduced" -p --factor=0.8
running parallel implementation
breadth-first traversal till k = 109
Liu Nguyen			[done]
../../bin/enumeration -n80 -s1.0 -q4093 -b2 -eln -i"reduced" -p --factor=0.8  3.00s user 0.01s system 100% cpu 2.997 total
```

## Test Solution ##
```
% diff solution_vector.dat samples_solution.dat
```
For the correct solution, this results in no difference.
