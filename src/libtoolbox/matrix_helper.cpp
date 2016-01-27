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
//     Implementation of some matrix helper operations.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>

#include <conversions.h>
#include <matrix_helper.h>

using namespace NTL;
using namespace std;

Mat<ZZ> RandomSqrMat(int n, int q)
{
	Mat<ZZ> m;
	m.SetDims(n, n);
	for (int i = 0; i < m.NumRows(); ++i)
	{
		m[i] = RandomVec(n, q);
	}
	return m;
}

/**
 * RandomVec
 * \brief returns a vector of length n drawn uniformly random from ZZ^n_q
 *
 * @params n
 * @params q
 */
Vec<ZZ> RandomVec(int n, long q)
{
	Vec<ZZ> v;
	ZZ mod = conv<ZZ>(q);

	v.SetLength(n);
	for (ZZ* it = v.begin(); it != v.end(); ++it)
	{
		RandomBnd(*it, mod);
	}
	return v;
}

/**
 * RemoveZeros
 * \brief remove zero vectors from matrix
 */
Mat<ZZ> RemoveZeros(Mat<ZZ> A)
{
	Mat<ZZ> newA;
	newA.SetDims(A.NumCols(), A.NumCols());
	int j = 0;
	for (int i = 0; i < A.NumRows(); i++)
	{
		if (!IsZero(A[i]))
		{
			newA[j] = A[i];
			j++;
		}
	}
	return newA;
}

Vec<RR> coeffs(Mat<ZZ> const& A, Vec<ZZ> const& t)
{
	RR det;
	Vec<RR> t_coeff;
	solve(det, t_coeff, conv<Mat<RR>>(A), conv<Vec<RR>>(t));

	return t_coeff;
}
vector<double> coeffs(matrix<long> const& A, vector<long> const& t)
{
    return to_stl<double>(coeffs(to_ntl<ZZ>(A), to_ntl<ZZ>(t)));
}

