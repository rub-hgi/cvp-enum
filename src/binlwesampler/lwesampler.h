// ----------------------------------------------------------------------------
// Title      : LWE Sampler
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : lwe_sampler.h
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

#ifndef __LWESAMPLER_H__
#define __LWESAMPLER_H__

#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>

#include <tuple>

std::tuple<NTL::Mat<NTL::ZZ>, NTL::Vec<NTL::ZZ>, NTL::Vec<NTL::ZZ>>
GenerateSamples(int n, long q, int m, double sigma, int tailfactor,
				int number_of_rects, int precision, int flag_binary,
				int flag_trinary, int flag_binary_secret, int flag_binary_a);

NTL::Vec<NTL::ZZ> RandomVec(int n, long q);
NTL::Mat<NTL::ZZ> PadMatrix(NTL::Mat<NTL::ZZ> A, long q, int m, int n);

NTL::Mat<NTL::ZZ> CreateLPerp(NTL::Mat<NTL::ZZ> const &A, int n, int m, long q);

#endif // __LWESAMPLER_H__
