// ----------------------------------------------------------------------------
// Title      : Enumeration Main
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : main_enumeration.cpp
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Main routine to controll the enumeration for closest lattice vector.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/ZZ.h>

#include <iostream>
#include <string>

#include <unistd.h>
#include <signal.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <conversions.h>
#include <io.h>
#include <matrix_helper.h>
#include <matrix_operations.h>
#include <signals.h>
#include <vector_templates.h>

#include "cmdline_enumeration.h"
#include "enumeration.h"

using namespace NTL;
using namespace std;

// command line arguments
double babaiBound_arg;
double delta_arg;
double factor_arg;
double factor_bin_arg;
long factor_lvl_arg;
double s_arg;
enum_dComp dComp_arg;
enum_rComp rComp_arg;
enum_enumeration enum_alg_arg;
int beta_arg;
int m_arg;
int n_arg;
int n_threads;
int parallel_flag;
int flag_binary_secret;
int flag_binary_a;
long q_arg;
string ifile_arg;
string ofile_arg;

bool got_sigterm = false;

// Internal Functions
void sigterm_handler(int _ignored);

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

	if (signal((int)SIGTERM, sigterm_handler) == SIG_ERR) {
		cerr << "failed to register SIGTERM handler, "
			 << "will not exit gracefully on SIGTERM" << endl;
	}
	got_sigterm = false;

	n_arg = args_info.dimension_arg;
	q_arg = (long)args_info.modulus_arg;
	s_arg = args_info.sigma_arg;
	beta_arg = args_info.beta_arg;
	enum_alg_arg = args_info.enumeration_arg;
	dComp_arg = args_info.dComp_arg;
	rComp_arg = args_info.rComp_arg;
	factor_arg = args_info.factor_arg;
	factor_bin_arg = args_info.factor_bin_arg;
	factor_lvl_arg = args_info.factor_lvl_arg;
	babaiBound_arg = args_info.babaiBound_arg;
	delta_arg = args_info.delta_arg;
	parallel_flag = args_info.parallel_flag;
	flag_binary_secret = args_info.binary_secret_flag;
	flag_binary_a = args_info.binary_a_flag;
	n_threads = args_info.n_threads_arg;
	ifile_arg = string(args_info.ifile_arg);
	ofile_arg = string(args_info.ofile_arg);

	cmdline_parser_free(&args_info);

	matrix<long> A = Read<matrix<long>>(ifile_arg + "_matrix.dat");
	const matrix<double> A_star = GSO(A);
	vector<long> v = Read<vector<long>>(ifile_arg + "_vector.dat");

	m_arg = (int)A.size();
	vector<long> d;
	vector<double> r;
	vector<long> t;

	if (parallel_flag) {
		switch (enum_alg_arg) {
		case enumeration_arg_lp:
		case enumeration_arg_ln: {
			cout << "running parallel implementation" << endl;
			break;
		}
		case enumeration__NULL:
		default:
			cerr << "requested parallel implementation not available for "
					"chosen enumeration algorithm" << endl;
			return EXIT_FAILURE;
		}
	}

	switch (enum_alg_arg) {
	case enumeration_arg_ntl: {
		t = NearestPlanesNTL(A, v);
		break;
	}

	case enumeration_arg_babai: {
		t = NearestPlanesBabaiOpt(A, v, A_star);
		break;
	}

	case enumeration_arg_lp: {
		switch (dComp_arg) {
		case dComp_arg_success: {
			d = ComputeD(VectorLengths(A_star), beta_arg, s_arg,
								 q_arg, factor_arg);
			break;
		}
		case dComp_arg_binary: {
			d = ComputeD_binary(VectorLengths(A_star), factor_arg, factor_bin_arg);
			break;
		}
		case dComp__NULL:
		default:
			cerr << "this should not happen" << endl;
			return EXIT_FAILURE;
		}

		cout << "d sequence used: " << endl << d << endl;
		if (parallel_flag) {
			size_t lvl = ComputeLvlNP(d, n_threads, factor_lvl_arg);

			t = NearestPlanesLPOptParall(A, A_star, d, v, q_arg, n_threads,
										 lvl);

		} else {
			t = NearestPlanesLPOpt(A, d, v, q_arg, A_star);
		}
		cout << "Lindner Peikert\t\t\t[done]" << endl;
		break;
	}

	case enumeration_arg_ln: {
		switch (rComp_arg) {
		case rComp_arg_length: {
			r = ComputeRlength(VectorLengths(A_star), s_arg, factor_arg, babaiBound_arg);
			break;
		}
		case rComp_arg_piece: {
			r = ComputeR_PiecewiseB(VectorLengths(A_star), s_arg, factor_arg);
			break;
		}
		case rComp__NULL:
		default:
			cerr << "this should not happen" << endl;
			return EXIT_FAILURE;
		}

		if (parallel_flag) {
			t = LengthPruningOptParall(A, A_star, r, v, q_arg, n_threads,
									   factor_lvl_arg);
		} else {
			t = LengthPruningOpt(A, A_star, r, v, q_arg);
		}
		cout << "Liu Nguyen\t\t\t[done]" << endl;
		break;
	}

	case enumeration__NULL:
	default:
		cerr << "Enumeration: this should not happen" << endl;
		return EXIT_FAILURE;
	}

	// reduce t modulo q, before writing it
	for (auto &i : t) {
		i = (i + q_arg) % q_arg;
	}

	if (flag_binary_a)
		Write(ofile_arg, (v - t) % q_arg);
	else if (flag_binary_secret) {
		t.resize(m_arg - n_arg);
		Write(ofile_arg, t);
	} else
		Write(ofile_arg, t);

	return EXIT_SUCCESS;
}

void sigterm_handler(int _ignored) {
	cout << "recieved SIGTERM " << _ignored
		 << ", stopping threads and write current best solution" << endl;
	got_sigterm = true;
}
