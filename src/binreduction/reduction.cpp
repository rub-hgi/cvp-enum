// ----------------------------------------------------------------------------
// Title      : Reduction
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : reduction.cpp
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Implementation of small wrapper around NTLs BKZ reduction
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#include <NTL/LLL.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>

#include <matrix_helper.h>

#include "reduction.h"

using namespace NTL;
using namespace std;

/**
 * ReduceMatrix
 * \brief reduces the matrix with BKZ and the given parameter
 *
 * @param
 */
Mat<ZZ> ReduceMatrix(Mat<ZZ> A, double delta, int beta, int prune) {
	LLLCheckFct CHECK = nullptr;
	const long VERBOSE = 1;
	// file to dump basis, 0 => no dump; initially 0
	LLLDumpFile = (char *)"BKZ_dump.dat";

	// seconds between status reports, initially 900s = 15min
	// dump every hour
	LLLStatusInterval = 3600;

	BKZ_RR(A, delta, beta, prune, CHECK, VERBOSE);

	// LLLDumpFile = nullptr;
	return RemoveZeros(A);
}
