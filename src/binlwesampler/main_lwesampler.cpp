// ----------------------------------------------------------------------------
// Title      : LWE Sampler Main
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : main_lwesampler.cpp
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Main routine to controll the LWE sampler from commandline, ties
//     everything together.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>
#include <NTL/LLL.h>

#include <iostream>
#include <tuple>
#include <string>

#include <io.h>
#include <matrix_operations.h>

#include "cmdline_lwesampler.h"
#include "lwesampler.h"

using namespace NTL;
using namespace std;

// command line arguments
int n_arg;
long q_arg;
int m_arg;
double s_arg;
string ofile_arg;
int flag_binary;
int flag_trinary;
int flag_binary_secret;
int flag_binary_a;

// sampling arguments - we do not change these anyway, so set them constant
const int TAILFACTOR = 13;
const int NUMBER_OF_RECTS = 63;
const int PRECISION = 106;

/**
 * main
 * \brief parses cli args, calls GenerateSamples and writes the result
 */
int main(int argc, char *argv[]) {
	gengetopt_args_info args_info;
	if (cmdline_parser(argc, argv, &args_info) != 0) {
		cerr << "failed parsing command line arguments" << endl;
		return EXIT_FAILURE;
	}

	n_arg = args_info.dimension_arg;
	q_arg = (long)args_info.modulus_arg;
	m_arg = args_info.samples_arg;
	s_arg = args_info.sigma_arg;
	ofile_arg = string(args_info.ofile_arg);
	flag_binary = args_info.binary_flag;
	flag_trinary = args_info.trinary_flag;
	flag_binary_secret = args_info.binary_secret_flag;
	flag_binary_a = args_info.binary_a_flag;

	tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> samples;
	// first generates LWE samples
	samples = GenerateSamples(n_arg, q_arg, m_arg, s_arg, TAILFACTOR,
							  NUMBER_OF_RECTS, PRECISION, flag_binary,
							  flag_trinary, flag_binary_secret, flag_binary_a);

	Mat<ZZ> A = get<0>(samples);

	if (flag_binary_a) {
		Vec<ZZ> w;
		A = transpose(A);
		LatticeSolve(w, A, get<1>(samples));
		Mat<ZZ> U;
		ZZ det;
		long r = image(det, A, U, 0);
		cout << "r = " << r << endl;
		Mat<ZZ> Ker = Kernel(U, r, q_arg);
		Write(ofile_arg + "_matrix.dat", Ker);

		for (int i = 0; i < m_arg; i++)
			w[i] = w[i] % q_arg;
		Write(ofile_arg + "_vector.dat", w);
	}
	// if the secret is binary, embed the target vector differently from normal
	// case
	else if (flag_binary_secret) {
		A = PadMatrix(get<0>(samples), q_arg, m_arg, n_arg);
		Mat<ZZ> LPerp = CreateLPerp(A, n_arg, m_arg, q_arg);
		Write(ofile_arg + "_matrix.dat", LPerp);

		Vec<ZZ> target = get<1>(samples);
		Vec<ZZ> ext_target;
		ext_target.SetLength(m_arg + n_arg);
		for (int i = 0; i < m_arg; i++)
			ext_target[i] = target[i];
		Write(ofile_arg + "_vector.dat", ext_target);
	} else {
		Write(ofile_arg + "_matrix.dat", PadMatrix(A, q_arg, m_arg, n_arg));
		Write(ofile_arg + "_vector.dat", get<1>(samples));
	}

	Write(ofile_arg + "_solution.dat", get<2>(samples));

	cmdline_parser_free(&args_info);
	return EXIT_SUCCESS;
}
