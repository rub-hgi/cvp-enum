// ----------------------------------------------------------------------------
// Title      : LWE Sampler
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : lwe_sampler.cpp
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Implementation of a lwesampler. It uses the ziggurat algorithm for
//     discrete gaussian sampling. The samples are returned as a matrix,
//     together with the target vector and the error free target.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>

#include <random>
#include <tuple>

#include <matrix_helper.h>

#include "lwesampler.h"
#include "ziggurat.h"

using namespace NTL;
using namespace std;

/**
 * GenerateSamples
 * \brief generates the requested number of samples and returns them in matrix
 *        form
 *
 * @param n length of each sample vector
 * @param q modulus
 * @param m number of samples
 * @param sigma for discrete gaussian
 * @param tailfactor for discrete gaussian
 * @param number_of_rects used in sampling algorithm (trade off speed/memory)
 * @param precision used in sampling algorithm
 */
tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> GenerateSamples(int n, long q, int m,
	double sigma, int tailfactor, int number_of_rects, int precision,
	int flag_binary, int flag_trinary, int flag_binary_secret, int flag_binary_a)
{
	Mat<ZZ> A;
	Vec<ZZ> t;
	Vec<ZZ> s;
	Vec<ZZ> e;
	Vec<ZZ> As; // solution vector -EK

	Vec<ZZ> u; //to test BinaryA
	Vec<ZZ> C1; //to test BinaryA
	Mat<ZZ> Abin;

	u.SetLength(m);
	C1.SetLength(n);
	Abin.SetDims(n, m);

	A.SetDims(m, n);  // m rows x n cols
	t.SetLength(m);
	As.SetLength(m);
	// generate LWE secret s
	if (flag_binary_secret)
	{
		s = RandomVec(n, 2);
	} else
	{
		s = RandomVec(n, q);
	}

	// generate error vector e
	if (flag_binary)
	{
		random_device rd;
		uniform_int_distribution<int> dist(0, 1);
		e.SetLength(m);
		for (ZZ* it = e.begin(); it != e.end(); ++it)
		{
			*it = ZZ(dist(rd));
		}
	// a uniform, trinary error (from {-1, 0, 1}) is used
	// in "Efficient Identity-Based Encryption over NTRU Lattices"
	// by Ducas, Lyubashevsky and Prest
	} else if (flag_trinary)
	{
		random_device rd;
		uniform_int_distribution<int> dist(-1, 1);
		e.SetLength(m);
		for (ZZ* it = e.begin(); it != e.end(); ++it)
		{
			*it = ZZ(dist(rd));
		}
	} else {
		Ziggurat sampler(number_of_rects, sigma, tailfactor, precision);
		e.SetLength(m);
		for (ZZ* it = e.begin(); it != e.end(); ++it)
		{
			*it = sampler.sample();
		}
	}

	// generate LWE samples
	tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> result;
	if (flag_binary_a)
	{
		u = RandomVec(m, 2);

		for (int i=0; i<n; i++)
		{
			Abin[i] = RandomVec(m, 2);
			InnerProduct(C1[i], Abin[i], u);
			C1[i] = C1[i] % q;
		}
		tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> result1(Abin, C1, u);
		result = result1;
	} else
	{
		for (int i=0; i<m; ++i)
		{
			A[i] = RandomVec(n, q);
			InnerProduct(t[i], A[i], s);
			As[i] = t[i] % q;
			t[i] += e[i];
			t[i] = t[i] % q;
		}
		tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> result1(A, t, As);
		result = result1;
	}

	//tuple<Mat<ZZ>, Vec<ZZ>, Vec<ZZ>> result(A, t, As);

	return result;
}

/**
 * PadMatrix
 * \brief padds an LWE-Matrix with qI matrix
 *
 * @params Mat<ZZ> - LWE matrix
 * @params q - LWE modulus
 */

Mat<ZZ> PadMatrix(Mat<ZZ> A, long q, int m, int n)
{
	Mat<ZZ> res;
	Mat<ZZ> qI;
	res.SetDims(m+n, m);

	A = transpose(A);

	for (int i = 0; i<n; i++)
	{
		res[i] = A[i];
	}

	diag(qI, conv<long>(m), conv<ZZ>(q));

	for (int i=n; i<m+n; i++)
	 {
		for (int j=0; j<m; j++)
		{
			res[i][j] = qI[i-n][j];
		}
	 }
	 return res;
}

Mat<ZZ> CreateLPerp(Mat<ZZ> const& A, int n, int m, long q)
{
	Mat<ZZ> LPerp;
	LPerp.SetDims(2*n+m, n+m);

	Mat<ZZ> qI;
	diag(qI, n+m, conv<ZZ>(q));

	for (size_t i=0; i<n; i++)
	{
		for (size_t j=0; j<n+m; j++)
		{
			if (j<m)
				LPerp[i][j] = -A[i][j] % q;
			if (j-i==m)
				LPerp[i][j] = 1;
		}
	}

	for (size_t i=n; i<2*n+m; i++)
	{
		for (size_t j = 0; j<n+m; j++)
			LPerp[i][j] = qI[i-n][j];
	}

	return LPerp;
}

