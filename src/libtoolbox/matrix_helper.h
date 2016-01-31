// ----------------------------------------------------------------------------
// Title      : Matrix Helper Functions
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : matrix_helper.h
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

#ifndef __MATRIX_HELPER_H__
#define __MATRIX_HELPER_H__

#include <NTL/matrix.h>
#include "NTL/mat_RR.h"
#include <NTL/vector.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <vector_templates.h>

NTL::Mat<NTL::ZZ> RandomSqrMat(int n, int q);
NTL::Vec<NTL::ZZ> RandomVec(int n, long q);
NTL::Mat<NTL::ZZ> RemoveZeros(NTL::Mat<NTL::ZZ> A);

NTL::Vec<NTL::RR> coeffs(NTL::Mat<NTL::ZZ> const &A,
						 NTL::Vec<NTL::ZZ> const &t);
std::vector<double> coeffs(matrix<long> const &A, std::vector<long> const &t);

#endif // __MATRIX_HELPER_H__
