// ----------------------------------------------------------------------------
// Title      : Conversions
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : conversions.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Templated implementations for conversions between NTL and stl.
//     This needs to be in a header file, as c++ does not support templated
//     implementations properly :/
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#ifndef __CONVERSIONS_H__
#define __CONVERSIONS_H__

#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <vector>

#include <vector_templates.h>

template <typename T, typename U>
matrix<T> to_stl(NTL::Mat<U> M)
{
	matrix<T> M_((size_t) M.NumRows(), std::vector<T>((size_t) M.NumCols()));

	for (int i = 0; i < M.NumRows(); ++i)
	{
		for (int j = 0; j < M.NumCols(); ++j)
		{
			M_[(size_t) i][(size_t) j] = NTL::conv<T>(M[i][j]);
		}
	}

	return M_;
}

template <typename T, typename U>
NTL::Mat<T> to_ntl(matrix<U> M)
{
	NTL::Mat<T> M_;
	if (M.size() == 0)
	{
		return M_;
	}

	M_.SetDims ((long) M.size(), (long) M[0].size());
	for (size_t i = 0; i < M.size(); ++i)
	{
		for (size_t j = 0; j < M[0].size(); ++j)
		{
			M_[(int) i][(int) j] = NTL::conv<T>(M[i][j]);
		}
	}
	return M_;
}

template <typename T, typename U>
std::vector<T> to_stl(NTL::Vec<U> u)
{
	std::vector<T> v((size_t) u.length());

	for (long i = 0; i < u.length(); ++i)
	{
		v[(size_t) i] = NTL::conv<T>(u[i]);
	}

	return v;
}

template <typename T, typename U>
NTL::Vec<T> to_ntl(std::vector<U> u)
{
	NTL::Vec<T> v;
	if (u.size() == 0)
	{
		return v;
	}

	v.SetLength ((long) u.size());
	for (size_t i = 0; i < u.size(); ++i)
	{
		v[(long) i] = NTL::conv<T>(u[i]);
	}
	return v;
}

#endif  // __CONVERSIONS_H__

