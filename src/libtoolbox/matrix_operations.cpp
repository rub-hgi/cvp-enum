// ----------------------------------------------------------------------------
// Title      : Reduction
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : reduction.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Implementation of some common matrix operations.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/matrix.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_RR.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

#include <cmath>
#include <numeric>
#include <vector>

#include <conversions.h>
#include <vector_templates.h>
#include <matrix_operations.h>

using namespace NTL;
using namespace std;

/**
 * GSO
 * \brief returns the Gram-Schmidt orthogonalized matrix
 *
 * @params B the basis
 */
matrix<double> GSO(matrix<long> const &B_std) {
	Mat<ZZ> B = to_ntl<ZZ>(B_std);
	Mat<RR> B_star;
	B_star.SetDims(B.NumRows(), B.NumCols());
	if (B.NumRows() > 0) {
		B_star[0] = conv<Vec<RR>>(B[0]);
		for (int i = 1; i < B.NumRows(); ++i) {
			B_star[i] = conv<Vec<RR>>(B[i]);
			for (int j = 0; j < i; ++j) {
				RR num; num = 0;
				RR denom; denom = 0;
				RR division; division = 0;

				InnerProduct(num, B_star[j], conv<Vec<RR>>(B[i]));
				InnerProduct(denom, B_star[j], B_star[j]);

				division = num / denom;

				B_star[i] = B_star[i] - division * B_star[j];
			}
		}
	}
	return to_stl<double>(B_star);
}

/**
 * GSO_norm
 * \brief returns the Gram-Schmidt orthonormalized matrix
 *
 * @params B the basis
 */
matrix<double> GSO_norm(matrix<long> const &B) {
	Mat<RR> B_star = to_ntl<RR>(GSO(B));

	RR inv_length_sq;
	for (int i = 0; i < B_star.NumRows(); ++i) {
		inv_length_sq = 1 / (B_star[i] * B_star[i]);
		B_star[i] = B_star[i] * inv_length_sq;
	}

	return to_stl<double>(B_star);
}

/**
 * muGSO
 * \brief on input B returnes a lower-triangular matrix mu such that B = mu*Q,
 *        with <b_i^*, b_i^*> / ||b_i^*|| = 1 on the main diagonal,
 *        Q - matrix with orthonormal rows.
 *
 * @params B the basis
 */
matrix<double> muGSO(matrix<long> const &B_std) {
	Mat<ZZ> B = to_ntl<ZZ>(B_std);
	Mat<RR> B_star;
	Mat<RR> muGSO;
	muGSO.SetDims(B.NumRows(), B.NumCols());
	B_star.SetDims(B.NumRows(), B.NumCols());
	if (B.NumRows() > 0) {
		B_star[0] = conv<Vec<RR>>(B[0]);
		InnerProduct(muGSO[0][0], B_star[0], B_star[0]);
		for (long i = 1; i < B.NumRows(); ++i) {
			B_star[i] = conv<Vec<RR>>(B[i]);
			for (long j = 0; j < i; ++j) {
				RR num; num = 0;
				RR denom; denom = 0;
				RR division; division = 0;

				InnerProduct(num, B_star[j], conv<Vec<RR>>(B[i]));
				InnerProduct(denom, B_star[j], B_star[j]);

				muGSO[i][j] = num / denom;
				B_star[i] = B_star[i] - muGSO[i][j] * B_star[j];
			}
			InnerProduct(muGSO[i][i], B_star[i], B_star[i]);
		}
	}
	return to_stl<double>(muGSO);
}

/**
 * muT
 * \brief on input GSO(B) and a vector t returnes the "embedded" vector t in the
 *        mu matrix
 *
 *        we use this trick to implement the optimized variant of NearestPlane
 *        as done in the NTL
 */
vector<double> muT(matrix<double> const &B_star_std,
				   vector<long> const &t_std) {
	Mat<RR> B_star = to_ntl<RR>(B_star_std);
	Vec<ZZ> t = to_ntl<ZZ>(t_std);
	long dim = t.length();
	Vec<RR> result;
	result.SetLength(dim);
	if (dim > 0) {
		for (long i = 0; i < dim; ++i) {
			RR num; num = 0;
			RR denom; denom = 0;

			InnerProduct(num, conv<Vec<RR>>(t), B_star[i]);
			InnerProduct(denom, B_star[i], B_star[i]);

			result[i] = num / denom;
		}
	}
	return to_stl<double>(result);
}

// to compute the coeff. of vector t with respect to basis B_star
vector<double> GSO_coeff(matrix<double> B_star, vector<long> t) {
	size_t m = t.size();

	Vec<RR> t_star;
	Vec<RR> t_conv = to_ntl<RR>(t);
	t_star.SetLength((int)m);

	mul(t_star, t_conv, to_ntl<RR>(B_star)); //--cases error TO DEBUG.

	//- to test wheter t_star*B_star = t --EK
	Vec<RR> tTest;
	tTest.SetLength((int)m);
	for (size_t i = 0; i < m; ++i) {
		tTest = tTest + t_star[(int)i] * to_ntl<RR>(B_star[i]);
	}

	return to_stl<double>(t_star);
}

/**
 * EuclideanNorm
 * \brief returns the euclidean / l2 Norm of its input vector
 *
 * @params u the vector
 */
double EuclideanNorm(vector<long> const &u) {
	double result = 0;

	for (auto const &i : u) {
		result += pow(abs(i), 2);
	}
	result = sqrt(result);

	return result;
}

/**
 * VectorLengths
 * \brief returns a vector with the lengths of each matrix vector
 *
 * @params M the matrix
 */
vector<double> VectorLengths(matrix<double> const &M) {
	vector<double> M_lengths;
	for (auto const &m_i : M) {
		M_lengths.push_back(
			inner_product(m_i.cbegin(), m_i.cend(), m_i.cbegin(), 0.0));
	}
	return M_lengths;
}

/**
 * Determinant
 * \brief wrapper for NTLs determinant implementation
 */
int Determinant(matrix<long> M) {
	return conv<int>(determinant(to_ntl<ZZ>(M)));
}

Mat<ZZ> Kernel(Mat<ZZ> U, long r, long q) {
	Mat<ZZ> ker;
	long m = U.NumCols();
	// ker.SetDims(m-r, m); //without q-ary part
	ker.SetDims(2 * m - r, m);
	for (int i = 0; i < m - r; i++)
		for (int j = 0; j < m; j++)
			ker[i][j] = U[i][j] % q;

	Mat<ZZ> qI;
	diag(qI, conv<long>(m), conv<ZZ>(q));

	for (int i = m - r; i < 2 * m - r; i++) {
		for (int j = 0; j < m; j++) {
			ker[i][j] = qI[i - m + r][j];
		}
	}

	return ker;
}
