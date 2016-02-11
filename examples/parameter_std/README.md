# Standard LWE Parameter #

## Generate Samples ##
```
% time ../../bin/lwesampler -n60 -m100 -s4 -q4093
../../bin/lwesampler -n60 -m100 -s4 -q4093  0.81s user 0.00s system 99% cpu 0.815 total
```

## Reduction ##
```
% time ../../bin/reduction -b5
---- LLL_RR status ----
elapsed time: 4:44, stage: 101, rank: 100, swaps: 204368
log of prod of lengths: 775.508
dumping to BKZ_dump.dat...
---- BKZ_RR status ----
elapsed time: 5:19, enum time: 0, iter: 5395
triv: 947, nontriv: 169, no ops: 4279, rank: 100, swaps: 214168
log of prod of lengths: 730.14
dumping to BKZ_dump.dat...
../../bin/reduction -b5  319.21s user 0.13s system 99% cpu 5:19.45 total

% cp samples_vector.dat reduced_vector.dat
```
Note that you have to copy the samples_vector file to reduced_vector file, because
the enumeration program looks for the matrix/vector in files with the same prefix.

## Enumeration ##
```
% time ../../bin/enumeration -n60 -s4 -q4093 -b5 -eln -i"reduced" -p --factor=1.1
running parallel implementation
breadth-first traversal till k = 93
Liu Nguyen			[done]
../../bin/enumeration -n60 -s4 -q4093 -b5 -eln -i"reduced" -p --factor=1.1  40.49s user 0.01s system 574% cpu 7.052 total
```

## Test Solution ##
```
% diff solution_vector.dat samples_solution.dat
```
For the correct solution, this results in no difference.
