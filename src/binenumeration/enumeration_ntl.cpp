#include <NTL/LLL.h>
#include <NTL/vec_ZZ.h>

#include <vector>

#include <conversions.h>

using namespace NTL;
using namespace std;

/**
 * NearestPlanesNTL
 * \brief takes a lattice basis and a target point and returns the 'relative
 *        closest' lattice point using NTL's implementation of the original
 *        (iterative) algorithm of Babai
 *
 * @params B the lattice basis
 * @params t the target point
 */
vector<long> NearestPlanesNTL(matrix<long> const &B, vector<long> const &t) {
	Vec<ZZ> v;
	NearVector(v, to_ntl<ZZ>(B), to_ntl<ZZ>(t));
	return to_stl<long>(v);
}
