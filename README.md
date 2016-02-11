# Parallel Implementation of BDD enumeration for LWE

## Compile and Install
The project depends on cmake (install cmake under Ubuntu) and
NTL. In order to compile and test this package, these dependencies are needed.

The following commands generate the build scripts, check for dependencies and
compile the source code.

```
	cd build
	cmake ..
	make
```

## Quickstart ##
After building the code, there will be four binaries in the bin-subfolder:
  1. bin/decoding
  2. bin/lwesampler
  3. bin/reduction
  4. bin/enumeration

The *decoding* binary runs a complete BDD enumeration attack on LWE, that is:
  * Generate the specified number of LWE samples
  * Reduce the resulting lattice basis
  * Run an enumeration on the reduced basis
  * Compare the resulting vector with the correct solution
  * Report runtimes

The three other binaries each do one step of the complete attack.
  * *lwesampler* generates an LWE instance and writes the matrix, target vector
    and solution to files.
  * *reduction* reads the matrix from a file and β-BKZ-reduces it.
  * *enumeration* reads the reduced matrix and vector from a file and runs an
    enumeration algorithm on them, writing the resulting solution vector to a
    file.

The solution can be checked by *diff*ing the outputted solution of lwesampler
and enumeration.

## Examples ##
The examples subdirectory contains example runs for each parameter set supported.
See the corresponding README.md for instructions on how to reproduce such runs.

## Structure of the Project ##

**TODO**:
  * Describe commandline parameter
