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
 *        and Peikert as proposed in "Elliptic Nearest Planes and the Complexity
 *        of Lattice-Based LWE Decoding"
 *
 * @params B the basis
 * @params s is used to blur the discrete structure of lattices over Z^n (?)
 * @params delta root hermite factor of reduction algorithm used before
 *         NearestPlanes
 * @params beta block size of BKZ run on basis
 */
vector<long> ComputeD(matrix<long> const &B, double s) {
	size_t m = B.size();
	vector<long> d(m);
	int det = Determinant(B);

	double mPow = pow(m, -1);
	double nom = pow(m, mPow);
	double denom = pow(abs(det), mPow * mPow);
	double delta = nom / denom;
	double y = (2 * s * sqrt(log2((double)m))) / ((double)m * sqrt(PI));

	for (size_t i = 0; i < m; ++i) {
		d[i] = (long)ceil(pow(delta, 2 * i) * y);
	}

	return d;
}
vector<long> ComputeD_success(matrix<double> A_mu, int beta, double s, int n_in,
							  long q, double factor) {
	size_t m = A_mu.size();
	vector<long> d(m);

	double tmp;

	for (size_t i = 0; i < m; ++i) {
		tmp = 3 * s * factor /
			  (2 * sqrt(A_mu[i][i])); // 3 might be reduced/changed --EK
		if (tmp < 1)
			d[i] = 1;
		else
			d[i] = (long)ceil(tmp);
	}
	return d;
}

/**
 * ComputeD_binary
 * compute the d-sequense for binary error
 */
vector<long> ComputeD_binary(matrix<double> A_mu, double s, int n,
							 double factor, double factor_bin) {
	size_t m = A_mu.size();
	vector<long> d(m);
	double tmp;

	// for the binary part
	for (size_t i = 0; i < m; i++) {
		if (i < n)
			tmp = factor_bin / (2 * sqrt(A_mu[i][i]));
		else
			tmp = 3 * s * factor / (2 * sqrt(A_mu[i][i]));
		if (tmp < 1)
			d[i] = 1;
		else
			d[i] = (long)ceil(tmp);
	}
	return d;
}

/**
 * ComputeR_LP
 * \brief TODO
 */
vector<double> ComputeR_LP(matrix<double> A_mu, vector<long> d) {
	size_t m = A_mu.size();
	vector<double> R(m + 1, 0);
	for (size_t i = 1; i <= m; ++i) {
		R[m - i] =
			R[m - i + 1] + (d[m - i]) * (d[m - i]) * A_mu[m - i][m - i] / 6;
		// cout << A_mu[i-1][i-1] << " ";
	}
	cout << endl;
	return R;
}

/**
 * ComputeRlength
 * \brief TODO
 */
vector<double> ComputeRlength(matrix<double> A_mu, double s, double factor,
							  double babaiBound) {
	size_t m = A_mu.size();
	vector<double> R(m);
	for (size_t i = 1; i <= m; ++i) {
		// babaiBound * ... is 'some' factor, maybe needs to be adjusted
		if (A_mu[m - i][m - i] > babaiBound * s * s) {
			// allow only one node on level i
			double prev_R_i = i == 1 ? 0 : R[m - i + 1];
			R[m - i] = prev_R_i + 0.25 * A_mu[m - i][m - i];
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
 *hardware
 */
size_t ComputeLvlNP(vector<long> const &d, int n_threads) {
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
		if (num_threads > number_of_max_threads)
			break;
	}

	// if num_threads == 1 then all d_i = 1 => no parallel runs
	if (num_threads == 1) {
		cout << "no parallel runs detected" << endl;
		lvl = 1;
	}

	cout << "number of maximal threads supported by hardware: "
		 << number_of_max_threads << endl
		 << "iteration level from where threads will be started: " << lvl
		 << endl << "number of threads to start at this level: " << num_threads
		 << endl << endl;
	return lvl;
}

/**
 * ComputeLvlLength
 * \brief computes the level from where to start parallel runs
 *
 * @params B the basis matrix
 * @params R
 * @params n_threads if we want to run less threads than available by the
 *hardware
 */
size_t ComputeLvlLength(matrix<long> const &B, vector<double> const &R,
						vector<long> t, int n_threads) {
	// check if we got a command line argument to bound the number of threads
	if ((n_threads > 0) && (n_threads <= number_of_max_threads)) {
		number_of_max_threads = n_threads;
	}

	// iterate down until number_of_max_threads is reached
	const matrix<double> B_star = GSO(B);
	const matrix<double> B_star_orth = GSO_norm(B);
	const vector<double> B_star_len = VectorLengths(B_star);
	size_t m = B.size();
	size_t lvl = 0;
	long num_threads = 1;
	double errL2 = 0.0;

	for (size_t i = R.size() - 1; i > 0; i--) {
		size_t level = R.size() - 1 - i;

		double interval =
			sqrt(fabs((R[m - level - 1] - errL2)) / B_star_len[m - level - 1]);
		double c_star = inner_product(t.cbegin(), t.cend(),
									  B_star_orth[m - level - 1].cbegin(), 0.0);
		long c = (long)ceil(c_star - interval);

		errL2 += pow(c_star - (double)c, 2) * B_star_len[m - level - 1];
		t -= c * B[m - level - 1];

		unsigned int nodes = ceil(interval);

		cout << "assume " << nodes << " nodes at level " << level << endl;
		if (nodes > 1) {
			num_threads *= nodes;
			lvl = level;
		}
		if (num_threads >= number_of_max_threads)
			break;
	}

	// if num_threads == 1 then all d_i = 1 => no parallel runs
	if ((num_threads == 1) || (lvl == R.size() - 1)) {
		cout << "no parallel runs detected" << endl;
		lvl = 1;
	}

	cout << "number of maximal threads supported by hardware: ";
	cout << number_of_max_threads << endl;
	cout << "iteration level from where threads will be started: ";
	cout << lvl << endl;
	cout << "number of threads to start at this level (assumend): ";
	cout << num_threads << endl << endl;

	return lvl;
}
