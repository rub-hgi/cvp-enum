# Examples of BDD Enumerations #

For each type of parameters (standard, binary and ternary) there is on
subfolder, including an example parameter set. These instances where created
on an *Intel Core i7-3820QM CPU @ 2.70GHz* and enumerations used all 8
available cores.

Please note that the runtimes vary. Also, using the standard success-factor
of 1.5 might lead to long runtimes, depending on the enumeration method used.
Thus, for tinkering with the implementation, it is best to start with a lower
factor, e.g. `--factor=0.8`, and increase it until the enumeration succeeds.

For detailed information on how the examples where generated, please see each
subfolders README.md.
