#include <iostream>
#include <vector>

#include <conversions.h>
#include <matrix_operations.h>
#include <vector_templates.h>

#include "enumeration.h"

using namespace std;

/**
 * NearestPlanesBabaiOpt
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice point using Babais original (iterative) algorithm
 *
 * @params B the lattice basis
 * @params t the target point
 */
vector<long> NearestPlanesBabaiOpt(matrix<long> const &B, vector<long> const &t,
								   matrix<double> const &B_star) {
	matrix<double> mu = muGSO(B);
	vector<double> mu_t = muT(B_star, t);
	vector<long> t_err(t);

	for (long i = (long)t.size() - 1; i >= 0; --i) {
		// c_star = round (muT_i)
		long c_star;
		if (mu_t[(size_t)i] >= 0)
			c_star = (long)ceil(mu_t[(size_t)i] - 0.5);
		else
			c_star = (long)floor(mu_t[(size_t)i] + 0.5);

		// update all muT_i's starting from the highest
		// muT_j = muT_j - c_star * muT_i
		for (long j = i - 1; j >= 0; --j)
			mu_t[(size_t)j] -= (double)c_star * mu[(size_t)i][(size_t)j];

		// row transforms
		t_err -= c_star * B[(size_t)i];
	}

	return t - t_err;
}
