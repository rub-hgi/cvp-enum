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
 * @params B_star is GSO(B)
 */
std::vector<long> NearestPlanesBabaiOpt(matrix<long> const &B,
										std::vector<long> const &t,
										matrix<double> const &B_star);

/**
 * NearestPlanesLP
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice points using Lindner Peikert's (LP's) algorithm
 *        http://www.cc.gatech.edu/~cpeikert/pubs/lwe-analysis.pdf
 *
 * @params B the (reduced) lattice basis
 * @params d determines the number of planes to project on
 * @params t the target point
 * @params q the modulus
 * @params B_star is GSO(B)
 */
std::vector<long> NearestPlanesLPOpt(matrix<long> const &B,
									 std::vector<long> const &d,
									 std::vector<long> const &t, long q,
									 matrix<double> const &B_star);
/**
 * NearestPlanesLPOptPaeallel
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice points using Lindner Peikert's (LP's) algorithm
 *        http://www.cc.gatech.edu/~cpeikert/pubs/lwe-analysis.pdf
 *        in parallel
 *
 * @params B the lattice basis
 * @params d determines the number of planes to project on
 * @params t the target point
 * @params q the modulus
 * @params B_star is GSO(B)
 * @params n_threads number of processors available
 * @params lvl level of the enumerarion tree to start parallel threads from.
 *             Determined by ComputeLvlNP
 */
std::vector<long> NearestPlanesLPOptParall(matrix<long> const &B,
										   matrix<double> const &B_star,
										   std::vector<long> const &d,
										   std::vector<long> t, long q,
										   int n_threads, size_t lvl);

/**
 * ComputeD
 * \brief returns a vector of d_i's for the NearestPlanes algorithm by Lindner
 *        and Peikert
 *
 * @params A_star_length lengths of vectors of Gram-Schmidt bassis
 * @params s Gaussian distribution parameter
 * @params q modulus
 * @params beta block size of BKZ run on basis
 * @params factor addition scaling factor to generate a sequence with
 *                larger/smaller success probability. By default, factor=1.0.
 */
std::vector<long> ComputeD(std::vector<double> A_star_length, int beta,
						   double s, long q, double factor);

/**
 * ComputeD
 * \brief returns a vector of d_i's for the NearestPlanes algorithm by Lindner
 *        and Peikert for binary/ternary error
 *
 * @params A_star_length lengths of vectors of Gram-Schmidt bassis
 * @params s Gaussian distribution parameter
 * @params q modulus
 * @params beta block size of BKZ run on basis
 * @params factor addition scaling factor to generate a sequence with
 *                larger/smaller success probability. By default, factor=1.0.
 */
std::vector<long> ComputeD_binary(std::vector<double> A_star_length,
								  double factor, double factor_bin);
std::vector<double> ComputeR_LP(matrix<double> A_star_length,
								std::vector<long> d);

/**
 * LengthPruningOpt
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice points using Pruned Enumeration algorithm
 *        https://hal.inria.fr/hal-00864361/PDF/LiuNguyen.pdf
 *
 * @params B the lattice basis
 * @params B_star = GSO(B)
 * @params R the bounds on squared error-lengths returned by ComputeRlength
 * @params t the targer
 * @params q the modulus
 */
std::vector<long> LengthPruningOpt(matrix<long> const &B,
								   matrix<double> const &B_star,
								   std::vector<double> const &R,
								   std::vector<long> const &t, long q);

/**
 * LengthPruning
 * \brief parallel version of LengthPruningOpt
 *
 * @params B the lattice basis
 * @params B_star = GSO(B)
 * @params R the bounds on squared error-lengths returned by ComputeRlength
 * @params t the targer
 * @params q the modulus
 * @params n_threads number of availabe processors
 * @params factor_lvl start the threads when nNodes(k) >= factor_lvl*n_threads,
 *                    where nNodes(k) is the number of nodes on level k
 */
std::vector<long> LengthPruningOptParall(matrix<long> const &B,
										 matrix<double> const &B_star,
										 std::vector<double> const &R,
										 std::vector<long> const &t, long q,
										 long nthreads, long factor_lvl);

/**
 * ComputeRlength
 * \brief computes the squared error-length allowed on each enumeration level.
 *        I.e. \|e[m-i] \|^2 <= R[m-i]. Note R[1] > R[2] > ... R[m]
 *        if the squared length of the Gram-Schmidt vectors > babai_bound * s^2,
 *        start Babai (one child allowed).
 *
 * @params A_star_length lengths of vectors of Gram-Schmidt bassis
 * @params s Gaussian distribution parameter
 * @params q modulus
 * @params factor additional scaling factor to generate a sequence with
 *                larger/smaller success probability. By default, factor=1.0.
 * @params babaiBound controls when to start Babai's CVP
 */
std::vector<double> ComputeRlength(matrix<double> A_star_length, double s,
								   double factor, double babaiBound);

/**
 * ComputeR_LP
 * \brief computes the squared error-length allowed on each enumeration level
 *        when the Lindner-Peikert enumeration is used. Currently, not in use.
 */
std::vector<double> ComputeR_PiecewiseB(std::vector<double> A_star_length,
										double s, double factor);

/**
 * ComputeLvl
 * \brief returns the lvl to start the parallel runs on
 *
 * @params d - the d-sequence returned by ComputeD
 * @params n_threads number of availabe processors
 * @params factor_lvl start the threads when nNodes(k) >= factor_lvl*n_threads,
 *                    where nNodes(k) is the number of nodes on level k
 */
size_t ComputeLvlNP(std::vector<long> const &d, int n_threads, long factor_lvl);

#endif // __ENUMERATION_H__
