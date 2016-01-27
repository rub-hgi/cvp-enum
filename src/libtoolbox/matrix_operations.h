// ----------------------------------------------------------------------------
// Title      : Matrix Operations
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : matrix_operations.h
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

#ifndef __MATRIX_OPERATIONS_H__
#define __MATRIX_OPERATIONS_H__

#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <vector>

template<typename T>
using matrix = std::vector<std::vector<T>>;

/**
 * GSO
 * \brief returns the Gram-Schmidt orthogonalized matrix
 *
 * @params B the basis
 */
matrix<double> GSO(matrix<long> const& B);

/**
 * GSO_norm
 * \brief returns the Gram-Schmidt orthonormalized matrix
 *
 * @params B the basis
 */
matrix<double> GSO_norm(matrix<long> const& B);

/**
 * muGSO
 * \brief on input B returnes a lower-trinagular matrix mu such that B = mu*Q,
 *        with ||b_i^*|| on the main diagonal, Q - matrix with orthonormal rows.
 *
 * @params B the basis
 */
matrix<double> muGSO (matrix<long> const& B);

std::vector<double> muT (matrix<double> const& B_star, std::vector<long> const& t);

//to compute the coeff. of vector t with respect to basis B_star
std::vector<double> GSO_coeff (matrix<double> B_star, std::vector<long> t);

/**
 * EuclideanNorm
 * \brief returns the euclidean / l2 Norm of its input vector
 *
 * @params u the vector
 */
double EuclideanNorm(std::vector<long> const& u);

/**
 * VectorLengths
 * \brief returns a vector with the lengths of each matrix vector
 *
 * @params M the matrix
 */
std::vector<double> VectorLengths (matrix<double> const& M);

/**
 * Determinant
 * \brief wrapper for NTLs determinant implementation
 */
int Determinant (matrix<long> M);

NTL::Mat<NTL::ZZ> Kernel (NTL::Mat<NTL::ZZ> U, long r, long q);

#endif  // __MATRIX_OPERATIONS_H__

