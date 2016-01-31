// ----------------------------------------------------------------------------
// Title      : Enumeration
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : enumeration.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Implementation of several enumeration algorithms.
//     NearestPlanesNTL is a simple wrapper around NTLs NearVector
//     implementation.
//     NearestPlanesBabaiOpt implements Babais NearVector algorithm.
//     NearestPlanesLPOpt implements Lindner Peikerts Nearest Planes algorithm.
//     PrunedEnumeration implements Liu Nguyen algorithm.
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#ifndef __ENUMERATION_H__
#define __ENUMERATION_H__

#include <NTL/ZZ.h>

/**
 * NearestPlanesNTL
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice point using NTL's implementation of the original
 *        (iterative) algorithm of Babai
 *
 * @params B the lattice basis
 * @params t the target point
 */
std::vector<long> NearestPlanesNTL(matrix<long> const &B,
								   std::vector<long> const &t);

/**
 * BabaiNearestPlanes
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice point using Babais original (iterative) algorithm
 *
 * @params B the lattice basis
 * @params t the target point
 */
std::vector<long> NearestPlanesBabaiOpt(matrix<long> const &B,
										std::vector<long> const &t,
										matrix<double> const &B_star);

/**
 * NearestPlanesLP
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice points using Lindner Peikert's (LP's) algorithm
 *        http://www.cc.gatech.edu/~cpeikert/pubs/lwe-analysis.pdf
 *        *Rec is a recursive, *Iter an iterative implementation.
 *
 * @params B the lattice basis
 * @params D precomputed gram-schmidt orthogonalization of B
 * @params d used to control recursion, in order to reduce error probability
 * @params t the target point
 * @params q the modulus
 */
std::vector<long> NearestPlanesLPOpt(matrix<long> const &B,
									 std::vector<long> const &d,
									 std::vector<long> const &t, long q,
									 matrix<double> const &B_star);

std::vector<long> NearestPlanesLPOptParall(matrix<long> const &B,
										   matrix<double> const &B_star,
										   std::vector<long> const &d,
										   std::vector<long> t, long q,
										   size_t lvl);

/**
 * PrunedEnumeration
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice points using Pruned Enumeration algorithm
 *        https://hal.inria.fr/hal-00864361/PDF/LiuNguyen.pdf
 *
 * @params B the lattice basis
 * @params B_star precomputed gram-schmidt orthogonalization of B
 * @params d used to control recursion, in order to reduce error probability
 * @params t the target point
 * @params q the modulus
 * @params lvl level from where to start parallel threads
 */
std::vector<long> PrunedEnumeration(matrix<long> const &B,
									matrix<double> const &B_star,
									std::vector<long> const &d,
									std::vector<long> t, long q, size_t lvl);

/**
 * ComputeD
 * \brief returns a vector of d_i's for the NearestPlanes algorithm by Lindner
 *        and Peikert as proposed in "Elliptic Nearest Planes and the Complexity
 *        of Lattice-Based LWE Decoding"
 *
 * @params B the basis
 * @params delta root hermite factor of reduction algorithm used before
 *               NearestPlanes
 * @params beta blocksize of BKZ reduction of B
 * @params s is used to blur the discrete structure of lattices over Z^n
 */
std::vector<long> ComputeD(matrix<long> const &B, double s);
std::vector<long> ComputeD_success(matrix<double> A_mu, int beta, double s,
								   int n_in, long q, double factor);
std::vector<long> ComputeD_binary(matrix<double> A_mu, double s, int n,
								  double factor, double factor_bin);
std::vector<double> ComputeR_LP(matrix<double> A_mu, std::vector<long> d);

/**
 * LengthPruning
 * \brief TODO
 *
 * @params B the lattice basis
 * @params t the target point
 * @params s the standard deviation of the sample-error
 * @params q the modulus
 */
std::vector<long> LengthPruningOpt(matrix<long> const &B,
								   std::vector<double> const &R,
								   std::vector<double> const &t, long q,
								   matrix<double> const &mu);
std::vector<long> LengthPruningOptParall(matrix<long> const &B,
										 matrix<double> const &B_star,
										 std::vector<double> const &R,
										 std::vector<long> const &t, long q,
										 size_t lvl);

/**
 * ComputeRlength
 * \brief returns a vector for computing the interval during length pruning
 *
 * @params m number of samples
 * @params s the standard deviation of the sample-error
 * @params factor to tune the number of iterations
 */
std::vector<double> ComputeRlength(matrix<double> A_mu, double s, double factor,
								   double babaiBound);

/**
 * ComputeLvl
 * \brief return the lvl from were to start the parallel runs
 */
size_t ComputeLvlNP(std::vector<long> const &d, int n_threads);
size_t ComputeLvlLength(matrix<long> const &B, std::vector<double> const &R,
						std::vector<long> t, int n_threads);

#endif // __ENUMERATION_H__
