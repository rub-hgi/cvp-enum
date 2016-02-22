#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

#include <conversions.h>
#include <matrix_operations.h>
#include <vector_templates.h>

#include "enumeration.h"

using namespace std;

static int number_of_max_threads = (int)thread::hardware_concurrency();

const double PI = 3.141592653589793238463;

/**
 * ComputeD
 * \brief returns a vector of d_i's for the NearestPlanes algorithm by Lindner
 *        and Peikert
 *
 * @params A_star_length lengths of vectors of Gram-Schmidt bassis
 * @params s Gaussian distribution parameter
 * @params q modulus
 * @params beta block size of BKZ run on basis
 * @params factor addition scaling factor to generate a sequence with
 *         larger/smaller success probability. By default, factor=1.0.
 */
vector<long> ComputeD(vector<double> A_star_length, int beta, double s, long q,
					  double factor) {
	size_t m = A_star_length.size();
	vector<long> d(m);

	double tmp;

	for (size_t i = 0; i < m; ++i) {
		tmp = 3 * s * factor /
			  (2 * sqrt(A_star_length[i])); // 3 might be reduced/changed --EK
		if (tmp < 1)
			d[i] = 1;
		else
			d[i] = (long)ceil(tmp);
	}
	return d;
}

/**
 * ComputeD_binary
 * \brief compute the d-sequense for binary error. Same as Compute_D above
 */
vector<long> ComputeD_binary(vector<double> A_star_length, double factor,
							 double factor_bin) {
	int m = (int)A_star_length.size();
	vector<long> d(m);
	double tmp;

	// for the binary part
	for (int i = 0; i < m; i++) {

		tmp = factor / (2 * sqrt(A_star_length[i]));
		if (tmp < 1)
			d[i] = 1;
		else
			d[i] = (long)ceil(tmp);
	}
	return d;
}

/**
 * ComputeRlength
 * \brief computes the squared error-length allowed on each enumeration level.
 *        I.e. \|e[m-i] \|^2 <= R[m-i]. Note R[1] > R[2] > ... R[m]
 *        if the squared length of the Gram-Schmidt vectors > babai_bound * s^2,
 *        start Babai (one child allowed).
 *
 * @params A_star_length lengths of vectors of Gram-Schmidt bassis
 * @params s Gaussian distribution parameter
 * @params q modulus
 * @params factor addition scaling factor to generate a sequence with
 *         larger/smaller success probability. By default, factor=1.0.
 * @params babaiBound controls when to start Babai's CVP
 */
vector<double> ComputeRlength(matrix<double> A_star_length, double s,
							  double factor, double babaiBound) {
	size_t m = A_star_length.size();
	vector<double> R(m);
	for (size_t i = 1; i <= m; ++i) {
		// babaiBound * ... is 'some' factor, maybe needs to be adjusted
		if (A_star_length[m - i][m - i] > babaiBound * s * s) {
			// allow only one node on level i
			double prev_R_i = i == 1 ? 0 : R[m - i + 1];
			R[m - i] = prev_R_i + 0.25 * A_star_length[m - i][m - i];
		} else {
			R[m - i] = (double)i * s * s * factor * factor;
		}
	}
	return R;
}

/**
 * ComputeLvlNP
 * \brief computes the level from where to start parallel runs - this is easy in
 *        the Lindner Peikert case, as the d vector codes the number of
 *        iterations at each level in the search tree. Thus we can just multiply
 *        d's elements, until we reach the number of desired threads
 *
 * @params d Lindner Peikerts d vector
 * @params n_threads if we want to run less threads than available by the
 *         hardware
 * @params factor_lvl is used to balance short running threads
 */
size_t ComputeLvlNP(vector<long> const &d, int n_threads, long factor_lvl) {
	// check if we got a command line argument to bound the number of threads
	if ((n_threads > 0) && (n_threads <= number_of_max_threads)) {
		number_of_max_threads = n_threads;
	}

	// compute which level in the iteration tree has enough paths for spawning
	// threads
	size_t lvl = 0;
	long num_threads = 1;
	for (size_t i = d.size() - 1; i > 0; i--) {
		if (d[i] > 1) {
			num_threads *= d[i];
			lvl = d.size() - 1 - i;
		}
		// factor_lvl balances short running threads
		if (num_threads > factor_lvl * number_of_max_threads)
			break;
	}

	// if num_threads == 1 then all d_i = 1 => no parallel runs
	if (num_threads == 1) {
		cout << "no parallel runs detected" << endl;
		lvl = 1;
	}

	cout << "number of maximal threads supported by hardware: "
		 << number_of_max_threads << endl;
	cout << "iteration level from where threads will be started: " << lvl
		 << endl;
	cout << "number of threads to start at this level: " << num_threads << endl
		 << endl;
	return lvl;
}

/**
 * ComputeRlength_Piesewise
 * \brief will be used to output the numerically optimized R-bounds for the
 *        length-pruning.
 */
vector<double> ComputeR_PiecewiseB(vector<double> A_star_length, double s,
								   double factor) {
	size_t m = A_star_length.size();
	vector<double> R_tmp(m + 1);
	vector<double> R(m);
	R_tmp[m] = 0;

	for (size_t i = 1; i <= m; ++i) {

		double a = s * s / A_star_length[m - i];
		cout << "a = " << a << " " << endl;
		if (a < 1)
			R[m - i] = (double)i * s * s * factor * factor;
		else if (a < 2)
			R[m - i] = (double)i * s * s * factor * factor * a;
		else
			R[m - i] = (double)i * s * s * factor * factor * a * .5;
	}

	// for (size_t i = 0; i < m; i++) {
	// 	R[i] = R_tmp[i];
	// }

	return R;
}
