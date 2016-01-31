// ----------------------------------------------------------------------------
// Title      : Reduction Main
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : main_reduction.cpp
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Main routine to controll the lattice reduciton.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/LLL.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>

#include <string>

#include <io.h>

#include "cmdline_reduction.h"
#include "reduction.h"

using namespace NTL;
using namespace std;

// command line arguments
int beta_arg;
double delta_arg;
string ifile_arg;
string ofile_arg;
int prune_arg;

/**
 * main
 * \brief parses cli args, calls ReduceMatrix and writes the result
 */
int main(int argc, char *argv[]) {
	gengetopt_args_info args_info;
	if (cmdline_parser(argc, argv, &args_info) != 0) {
		cerr << "failed parsing command line arguments" << endl;
		return EXIT_FAILURE;
	}

	beta_arg = args_info.beta_arg;
	delta_arg = args_info.delta_arg;
	ifile_arg = string(args_info.ifile_arg);
	ofile_arg = string(args_info.ofile_arg);
	prune_arg = args_info.prune_arg;

	Mat<ZZ> A = Read<Mat<ZZ>>(ifile_arg);
	Mat<ZZ> A_red = ReduceMatrix(A, delta_arg, beta_arg, prune_arg);

	Write(ofile_arg, A_red);

	cmdline_parser_free(&args_info);
	return EXIT_SUCCESS;
}
