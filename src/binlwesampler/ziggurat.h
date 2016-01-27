// ----------------------------------------------------------------------------
// Title      : Ziggurat Sampler
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : ziggurat.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     This ziggurat implementation is an only sligthly modified version of
//     https://eprint.iacr.org/2013/510.pdf - the sampler is wrapped by a
//     c++ class. The original implementation is available at
//     https://www.cdc.informatik.tu-darmstadt.de/~pschmidt/
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#ifndef __ZIGGURAT_H__
#define __ZIGGURAT_H__

#include <NTL/RR.h>
#include <NTL/ZZ.h>

class Ziggurat {
	public:
		Ziggurat(int, double, int, int);
		~Ziggurat();
		NTL::ZZ sample();
	private:
		int omega; // precision
		NTL::RR m;
		NTL::RR sigma;
		NTL::RR tailcut;
		NTL::RR *rectys;
		NTL::ZZ *rectxs;
};

#endif  // __ZIGGURAT_H__
