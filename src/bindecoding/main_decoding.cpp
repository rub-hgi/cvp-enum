// ----------------------------------------------------------------------------
// Title      : Decoding Main
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : main_decoding.cpp
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Main routine to controll the complete decoding process.
//     This includes every step, from generating LWE samples, reducing the
//     resulting matrix and finally enumerate for the closest lattice vector.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/matrix.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>

#include <cassert>

#include <chrono>
#include <iostream>
#include <map>
#include <numeric>
#include <tuple>

#include <unistd.h>
#include <signal.h>

#include <signals.h>
#include <conversions.h>
#include <io.h>
#include <matrix_helper.h>
#include <matrix_operations.h>
#include <vector_templates.h>

#include "enumeration.h"
#include "lwesampler.h"
#include "reduction.h"
#include "cmdline_decoding.h"

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
int flag_binary;
int flag_trinary;
int flag_binary_secret;
int flag_binary_lwe;
int flag_binary_sis;
long q_arg;
int prune_arg;

// sampling arguments - we do not change these anyway, so set them constant
const int TAILFACTOR = 13;
const int NUMBER_OF_RECTS = 63;
const int PRECISION = 106;

bool got_sigterm = false;

// Internal Functions
void sigterm_handler(int _ignored);

// Internal Functions
tuple<matrix<long>, vector<long>, vector<long>>
Sample(chrono::duration<double> &duration);
matrix<long> Reduce(matrix<long> const &A, chrono::duration<double> &duration);
vector<long> Enumerate(matrix<long> const &A, vector<long> const &v,
					   chrono::duration<double> &duration);
bool Check(vector<long> const &solution, vector<long> const &t);

/**
 * main
 * \brief parses cli args, generates LWE samples and run BDD attack on it
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
	m_arg = args_info.samples_arg;
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
	n_threads = args_info.n_threads_arg;
	flag_binary = args_info.binary_flag;
	flag_trinary = args_info.trinary_flag;
	flag_binary_secret = args_info.binary_secret_flag;
	flag_binary_lwe = args_info.binary_lwe_flag;
	flag_binary_sis = args_info.binary_sis_flag;
	prune_arg = args_info.prune_arg;

	cmdline_parser_free(&args_info);

	// timing
	chrono::duration<double> duration_sample;

	// sampling
	tuple<matrix<long>, vector<long>, vector<long>> samples =
		Sample(duration_sample);

	vector<long> solution = get<2>(samples);
	vector<long> t;

	chrono::duration<double> duration_reduction;
	chrono::duration<double> duration_enumeration;

	// reduction
	matrix<long> A_red = Reduce(get<0>(samples), duration_reduction);

	// enumeration
	t = Enumerate(A_red, get<1>(samples), duration_enumeration);

	cout << endl << "time for sampling:\t\t" << duration_sample.count() << " s"
		 << endl << "time for reduction:\t\t" << duration_reduction.count()
		 << " s" << endl << "time for enumeration:\t\t"
		 << duration_enumeration.count() << " s" << endl
		 << "\ntime all together:\t\t"
		 << (duration_sample + duration_reduction + duration_enumeration)
				.count() << " s" << endl;

	cout << "error vector = " << (get<1>(samples) - t) % q_arg << endl;

	if (!Check(solution, t)) {
		cerr << "error!" << endl << "solution = " << solution << endl
			 << "found t  = " << t << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

// Internal Function Implementation

/**
 * Sample
 * \brief
 */
tuple<matrix<long>, vector<long>, vector<long>>
Sample(chrono::duration<double> &duration) {

	tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> samples;

	chrono::time_point<chrono::system_clock> start;
	start = chrono::system_clock::now();

	samples =
		GenerateSamples(n_arg, q_arg, m_arg, s_arg, TAILFACTOR, NUMBER_OF_RECTS,
						PRECISION, flag_binary, flag_trinary,
						flag_binary_secret, flag_binary_lwe, flag_binary_sis);
	// matrix<long> A = to_stl<long>(get<0>(samples));
	Mat<ZZ> A = get<0>(samples);
	vector<long> v;
	vector<long> solution = to_stl<long>(get<2>(samples));

	//-------For Binary-SIS---------
	if (flag_binary_sis) {
		Vec<ZZ> w;
		A = transpose(A);
		LatticeSolve(w, A, get<1>(samples));
		Mat<ZZ> U;
		ZZ det;
		long r = image(det, A, U, 0);
		cout << "r = " << r << endl;
		A = Kernel(U, r, q_arg);

		v = to_stl<long>(w) % q_arg;
	}
	// if the secret is binary, embed the target vector differently from normal
	// case
	else if (flag_binary_secret == 1) {
		cout << "the secret is binary" << endl;

		A = PadMatrix(A, q_arg, m_arg, n_arg);
		A = CreateLPerp(A, n_arg, m_arg, q_arg, s_arg);

		v.resize(m_arg + n_arg);
		for (int i = 0; i < m_arg; i++)
			v[i] = conv<long>(get<1>(samples)[i]);
	} else {
		A = PadMatrix(A, q_arg, m_arg, n_arg);
		v = to_stl<long>(get<1>(samples));
	}

	chrono::time_point<chrono::system_clock> end;
	end = chrono::system_clock::now();

	duration = end - start;
	return tuple<matrix<long>, vector<long>, vector<long>>(to_stl<long>(A), v,
														   solution);
}

/**
 * Reduce
 * \brief TODO
 */
matrix<long> Reduce(matrix<long> const &A, chrono::duration<double> &duration) {
	chrono::time_point<chrono::system_clock> start;
	chrono::time_point<chrono::system_clock> end;

	matrix<long> A_red;

	start = chrono::system_clock::now();
	A_red = to_stl<long>(
		ReduceMatrix(to_ntl<ZZ>(A), delta_arg, beta_arg, prune_arg));
	end = chrono::system_clock::now();

	cout << "Reduction\t\t\t[done]" << endl;
	duration = end - start;
	return A_red;
}

/**
 * Enumerate
 * \brief TODO
 */
vector<long> Enumerate(matrix<long> const &A, vector<long> const &v,
					   chrono::duration<double> &duration) {
	chrono::time_point<chrono::system_clock> start;
	chrono::time_point<chrono::system_clock> end;

	const matrix<double> A_star_orth = GSO_norm(A);
	const matrix<double> A_star = GSO(A);
	const matrix<double> A_mu = muGSO(A);

	m_arg = (int)A.size();
	vector<long> d;
	vector<double> r;
	vector<long> t;

	switch (enum_alg_arg) {
	case enumeration_arg_ntl: {
		start = chrono::system_clock::now();
		t = NearestPlanesNTL(A, v);
		end = chrono::system_clock::now();

		cout << "Babai's NearestPlane (NTL)\t[done]" << endl;
		break;
	}

	case enumeration_arg_babai: {

		start = chrono::system_clock::now();
		t = NearestPlanesBabaiOpt(A, v, A_star);
		end = chrono::system_clock::now();

		cout << "Babai's NearestPlane (Opt)\t[done]" << endl;
		break;
	}

	case enumeration_arg_lp: {
		switch (dComp_arg) {
		case dComp_arg_success: {
			d = ComputeD(VectorLengths(A_star), beta_arg, s_arg, q_arg,
						 factor_arg);
			break;
		}

		case dComp_arg_binary: {
			d = ComputeD_binary(VectorLengths(A_star), factor_arg,
								factor_bin_arg);
			break;
		}
		case dComp__NULL:
		default:
			cerr << "this should not happen" << endl;
			return t;
		}

		cout << "d sequence used: " << endl << d << endl;
		if (parallel_flag) {
			size_t lvl = ComputeLvlNP(d, n_threads, factor_lvl_arg);

			start = chrono::system_clock::now();
			t = NearestPlanesLPOptParall(A, A_star, d, v, q_arg, n_threads,
										 lvl);
			end = chrono::system_clock::now();
		} else {
			start = chrono::system_clock::now();
			t = NearestPlanesLPOpt(A, d, v, q_arg, A_star);
			end = chrono::system_clock::now();
		}
		cout << "Lindner Peikert\t\t\t[done]" << endl;
		break;
	}

	case enumeration_arg_ln: {
		switch (rComp_arg) {
		case rComp_arg_length: {
			r = ComputeRlength(A_mu, s_arg, factor_arg, babaiBound_arg);
			break;
		}
		case rComp_arg_piece: {
			r = ComputeR_PiecewiseB(VectorLengths(A_star), s_arg, factor_arg);
			break;
		}
		case rComp__NULL:
		default:
			cerr << "this should not happen" << endl;
			return t;
		}

		if (parallel_flag) {
			start = chrono::system_clock::now();
			t = LengthPruningOptParall(A, A_star, r, v, q_arg, n_threads,
									   factor_lvl_arg);
			end = chrono::system_clock::now();

		} else {
			start = chrono::system_clock::now();
			t = LengthPruningOpt(A, A_star, r, v, q_arg);
			end = chrono::system_clock::now();
		}
		cout << "Liu Nguyen\t\t\t[done]" << endl;
		break;
	}

	case enumeration__NULL:
	default:
		cerr << "Enumeration: this should not happen" << endl;
		return t;
	}

	duration = end - start;
	return t % q_arg;
}

void sigterm_handler(int _ignored) {
	cout << "recieved SIGTERM" << _ignored
		 << ", stopping threads and write current best solution" << endl;
	got_sigterm = true;
}

/**
 * Check
 * \brief checks the first size(solution) entries: solution == t
 */
bool Check(vector<long> const &solution, vector<long> const &t) {
	for (size_t i = 0; i < solution.size(); ++i)
		if (solution[i] != t[i]) {
			cout << "solution[" << i << "] = " << solution[i] << " != " << t[i]
				 << " = t[" << i << "]" << endl;
			return false;
		}
	return true;
}
