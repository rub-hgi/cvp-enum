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

The project is split in four executables and one convenience library. Each part
is described in the following sections.

### libtoolbox ###
Convenience library, covering matrix, io, conversion functions, which are needed
in each of the four binaries.

### lwesampler ###
Generates an LWE instance. Uses the discrete Ziggurat implementation from
[here](https://eprint.iacr.org/2013/510), the code is available
[here](https://www.cdc.informatik.tu-darmstadt.de/~pschmidt/implementations/ziggurat/ziggurat-src.zip).
Several options are available to customize the generated LWE instance:

```
% bin/lwesampler -h
lwesampler 0.0.1

The program generates the specified number of LWE samples.

Usage: lwesampler [OPTIONS]...

  -h, --help            Print help and exit
  -V, --version         Print version and exit
  -n, --dimension=INT   lattice dimension for LWE samples
  -q, --modulus=LONG    lattice modulus
  -m, --samples=INT     number of samples
  -s, --sigma=DOUBLE    standard deviaton of error distribution


Optional Arguments:
  the following options can be used to tune internal behaviour of program
  -o, --ofile=FILENAME  output file  (default=`samples')
  -2, --binary          generate error from binary  uniform distribution
                          {0, 1}  (default=off)
  -3, --trinary         generate error from trinary uniform distribution {-1,
                          0, 1}  (default=off)
      --binary_secret   generate secrect from binary uniform distribution {0,
                          1}  (default=off)
      --binary_lwe      generate A from binary uniform distribution {0, 1}
                          (default=off)
      --binary_sis      generate A from binary uniform distribution {0, 1}
                          (default=off)
```

All of these options should be self-explaining.

### reduction ###
Mainly a wrapper around NTL's BKZ implementation.

```
% bin/reduction -h
reduction 0.0.1

Reads a basis from a file and run BKZ reduction with given parameters.

Usage: reduction [OPTIONS]...

  -h, --help            Print help and exit
  -V, --version         Print version and exit
  -b, --beta=INT        blocksize of BKZ reduction


Optional Arguments:
  the following options can be used to tune internal behaviour of program
  -d, --delta=DOUBLE    delta of BKZ reduction  (default=`0.99')
  -i, --ifile=FILENAME  input file  (default=`samples_matrix.dat')
  -o, --ofile=FILENAME  output file  (default=`reduced_matrix.dat')
  -P, --prune=INT       pruning factor for BKZ reduction  (default=`0')
```

Again, the arguments should be self-explaining. ```beta```, ```delta```, and ```prune```
are used to tune the behaviour of the BKZ implementation. If the input/output files have
non-standard names, ```ifile``` and ```ofile``` needs to be adjusted.

### enumeration ###
Runs the choosen enumeration algorithm.

TODO: describe enumeration implementations

```
% bin/enumeration -h
enumeration 0.0.1

Reads a basis from a file and a target vector, runs an enumeration algorithm
on it and outputs the result.

Usage: enumeration [OPTIONS]...

  -h, --help                   Print help and exit
  -V, --version                Print version and exit
  -n, --dimension=INT          lattice dimension for LWE samples
  -q, --modulus=LONG           lattice modulus
  -s, --sigma=DOUBLE           standard deviaton of error distribution
  -b, --beta=INT               blocksize of BKZ reduction
  -e, --enumeration=ENUM       which enumeration algorithm to use (lp = Lindner
                                 Peikert, ln = Liu Nguyen)  (possible
                                 values="ntl", "babai", "lp", "ln"
                                 default=`ln')


Optional Arguments:
  the following options can be used to tune internal behaviour of program
      --babaiBound=DOUBLE      when running LengthPruning, this factor
                                 controlls the generation of the R sequecne and
                                 especially when to switch to Babais
                                 enumeration  (default=`4')
      --dComp=ENUM             how to compute d Sequence for LP's enumeration
                                 (possible values="success", "binary"
                                 default=`success')
      --rComp=ENUM             how to compute R Sequence for Length Pruning
                                 (possible values="length", "piece"
                                 default=`length')
      --factor=DOUBLE          controls the number of iterations done during
                                 enumeration  (default=`1.5')
      --factor_bin=DOUBLE      controls the number of iterations done when
                                 using binary secret d sequences
                                 (default=`1.0')
      --factor_lvl=LONG        used to balance short running threads by
                                 increasing the overall number of threads
                                 (default=`8')
  -d, --delta=DOUBLE           delta of BKZ reduction  (default=`0.99')
  -p, --parallel               run parallel implementations of enumeration
                                 algorithms  (default=off)
      --n-threads=INT          number of threads to use in parallel
                                 implementations  (default=`0')
  -i, --ifile=FILENAME PREFIX  input file prefix, "_{matrix,vector}.dat" will
                                 be appended.  (default=`samples')
  -o, --ofile=FILENAME         output file  (default=`solution_vector.dat')
      --binary_secret          generate secrect from binary uniform
                                 distribution {0, 1}  (default=off)
      --binary_a               generate A from binary uniform distribution {0,
                                 1}  (default=off)
```

```enumeration``` argument chooses the used algorithm:
  * **NTL** is the Babai implementation from the NTL library (called NearVector)
  * **lp** is Lindner-Peikerts Nearest Planes
  * **ln** is Liu-Nquyens Linear Length Pruning

Several other arguments can be used to fine-tune the algorithms success-probability
and runtime. ```dComp``` and ```rComp``` enables different functions to compute the
d sequence (Lindner Peikert) or r sequence (Liu Nguyen). The ```factor``` argument
controls how big these sequences get and thus controls the success-probability
(higher factor, higher success-probabilitiy and higher runtime). For the length
pruning ```factor=0.8``` seems to be a good starting point, Lindner-Peikert can be
a bit higher (one might like to check the d-sequence outputted for Lindner-Peikert).

When one notices less CPU usagen than expectet (say ```n-threads=8``` but only four
cores do work), short running threads might be the reason. Thus some threads finished
early and there are no more threads to start. ```factor_lvl``` can be used to start
spawning threads from a lower iteration level, where more threads to start are
available.

Note that parallel runs are switched off by default (```parallel``` flag).
