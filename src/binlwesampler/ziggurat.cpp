// ----------------------------------------------------------------------------
// Title      : Ziggurat Sampler
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : ziggurat.cpp
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

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <string.h>

#include "ziggurat.h"

using namespace std;
using namespace NTL;

//#define PRINT

////////////////////////////////////////////////////////////////////////////////
// internal functions //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int rejSampleBits(int);
void rejSample(ZZ*, ZZ);
void rejSample(int*, int);
void rejSampleRect(int*, int);
void rho(RR*, RR*, RR*);
RR rho(RR, RR);
RR compute_recurrence(RR*, RR*, RR, RR, RR, int);

////////////////////////////////////////////////////////////////////////////////
// implementation of external functions ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// - m = number of rectangles
// - sigma = "std. deviation" of (discrete Gaussian) distribution
// - precision of the computation
Ziggurat::Ziggurat(int number_of_rects, double s, int tailfactor, int precision) :
    omega(precision),
    m(ConvPrec(number_of_rects, omega)),
    sigma(ConvPrec(s, precision)),
    tailcut(to_RR(sigma * tailfactor)),
    rectys(new RR[number_of_rects+1]), rectxs(new ZZ[number_of_rects+1])
{
    RR y0;

    // set Seed for NTL-internal PRNG
    fstream urandom("/dev/urandom", ios::in|ios::binary);
    unsigned char randarray[512];
    urandom.read((char*)randarray,512);
    urandom.close();
    SetSeed(ZZFromBytes(randarray,512));

    // set precisions
    // max. approximation error
    RR prec = power2_RR(-omega);
    RR prec2 = prec * prec;
    // output-precision of numbers
    RR::SetOutputPrecision(omega);

    // initialize different variables
    RR bestdiff = to_RR(3);
    RR *xis = new RR[number_of_rects+1];
    RR *bestxis = new RR[number_of_rects+1];
    for (int i = 0; i <= number_of_rects; ++i) {
        bestxis[i] = -1;
    }
    RR *rhos = new RR[number_of_rects+1];
    RR c = 1 + 1 / m;

    if (m == 1) {
        rectxs[0] = to_ZZ(1);
        rectys[0] = rho(to_RR(1), sigma);
        rectxs[1] = to_ZZ(0);
        rectys[1] = rho(to_RR(0), sigma);

        y0 = to_RR(1);
        rectys[0] = y0;
    } else {
        // m != 1
        // increase right bound until reaching 14*sigma and compute possible partition(s)
        while (tailcut < to_RR((tailfactor + 1) * sigma))
        {
            xis[number_of_rects] = tailcut;
            RR cu = to_RR(0);
            RR cl = to_RR(1);
            RR cc;
            RR difference = to_RR(-1);
            RR lastdiff = to_RR(-2);

            // try to minimize the distance of y0 to 1 (y0 = y-value to x=0)
            while (difference < 0 ||
                    (abs(difference) > prec && abs(difference-lastdiff) > prec)
                  ) {
                cc = c;
                lastdiff = difference;
                difference = compute_recurrence(xis, rhos, c, m, sigma, tailfactor) - to_RR(1);
                if (difference == -2) {
                    break; // in case of any failure in partition-computation
                }
                if (difference >= 0) {
                    // if possible partition found, renew best solution,
                    // since difference to 1 is smaller than before
                    for (int i = 0; i <= number_of_rects; ++i) {
                        bestxis[i] = xis[i];
                    }
                    cc = c;
                    cu = c;
                    bestdiff = difference;
                } else {
                    cl = c;
                }
                // do some tricks with the "area" c in order to improve solution
                if (cu < cl) {
                    c += to_RR(1) / m;
                } else {
                    c = (cu + cl) /to_RR(2);
                }
                if (c >= 11) {
                    break;
                }
                if (difference == lastdiff) {
                    break;
                }
            }
            // if while-loop did not improve anything, increase tailcut and go to next iteration
            if (difference < 0 ||
                (abs(difference) > prec && abs(difference-lastdiff) > prec)
                ) {
                tailcut++;
            } else {
            // else stop the partition-computation
                break;
            }
        }
    }

    // if valid solution exists, output it to terminal and return 0 (success)
    if (bestxis[number_of_rects] != -1) {
        // output precision, number of rectangles m, (std. deviation) sigma, and y0>=1 (the value for x=0)
        y0 = to_RR(1) + bestdiff;

        // output the "x_i's", i.e. values on the x-axis
        for (int i = 0; i <= number_of_rects; ++i) {
            rectxs[i] = TruncToZZ(bestxis[i]);
            rectys[i] = rho(bestxis[i], sigma);
        }
        rectys[0] = y0;
    }
}

Ziggurat::~Ziggurat()
{
    delete[] rectxs;
    delete[] rectys;
}

ZZ Ziggurat::sample() {
    ZZ x;

    long outPrecision = 30;

    // set precisions
    RR::SetPrecision(omega); // max. approximation error
    RR::SetOutputPrecision(outPrecision); // output-precision of numbers

    int i; // rectangle i <- {1,...,m}
    int s = 1 - 2 * rejSampleBits(1); // sign <- +/- 1
    ZZ xurb;
    //~ int debug_number_rejections = 0, debug_right_rect = 0;
    while (true) {
        rejSampleRect(&i, to_int(m));
        i++; // sample rectangle uniformly
        xurb = rectxs[i] + 1;
        rejSample(&x, xurb); // res <- [0, \floor{x_{i}}]
        if (x != 0 && x <= rectxs[i-1]) {
            break;
        } else {
            // the case x=0 is special due to 0=-0, therefore we have to
            // halve the prob. for 0 which results in 1/2; so with p=1/2
            // p gets accepted and with p=1/2 rejected -> "coin toss"
            if (x == 0) {
                if (rejSampleBits(1) == 0) {
                    break;
                }
            } else {
                ZZ y;
                ZZ yprime;
                ZZ yfactor = TruncToZZ(power2_RR(omega));
                rejSample(&yprime, yfactor);
                y = TruncToZZ(to_RR(yprime) * (rectys[i-1] - rectys[i]));
                    // y <- [0,\floor{2^{\omega} (\rho_{\sigma}(x_{i-1}) - \rho_{\sigma}(x_{i}))}]
                if (y <= TruncToZZ(to_RR(yfactor) * (rho(to_RR(x), sigma) - rectys[i]))) {
                    break;
                }
            }
        }
    }
    x *= s;

    return x; //+ to_ZZ(tailcut);
}

////////////////////////////////////////////////////////////////////////////////
// implementation of internal functions ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int rejSampleBits(int m) {return (int) RandomBits_long(m);}
void rejSample(ZZ* res, ZZ m) {*res = RandomBnd(m);}
void rejSample(int* res, int m) {*res = (int) RandomBnd(m);}
void rejSampleRect(int* res, int m) {*res = (int) RandomBnd(m);}

/** computes e^{-1/2 (x/sigma)^2} **/
void rho(RR* res, RR* x, RR* sigma) {*res = exp(-power(*x / *sigma, 2)/2);}
RR rho(RR x, RR sigma) {return exp(-power(x / sigma, 2)/2.0);}

//- xis = storage for x-values
//- rhos = storage for y-values
//- c = area of rectangles
//- sigma = "std. deviation" of (discrete Gaussian) distribution
//- m = number of rectangles
RR compute_recurrence(RR* xis, RR* rhos, RR c, RR m, RR sigma, int tailfactor) {
    // m = #rectangles, thus we need m+1 x-values ("x_i's" or "xis")
    int em = to_int(m) + 1;

    // helping variables for computation
    RR om = to_RR(1) / to_RR(m);
    RR m2 = to_RR(-2);
    RR o2 = to_RR(1) / to_RR(2);

    // area of each rectangle
    RR area = sigma * om * sqrt(ComputePi_RR() / to_RR(2)) * to_RR(c);

    // set the largest x-value to 13*sigma
    xis[em-1] = to_RR(tailfactor * sigma);

    // rhos = corresponding y-values to the x_i's
    // manually set y-value for largest x-value (rounded to the next-largest integer)
    rhos[em-1] = rho(to_RR(TruncToZZ(xis[em-1]) + to_ZZ(1)), sigma);
    RR sqrtv = m2 * log(area/to_RR(TruncToZZ(xis[em-1]) + to_ZZ(1)));
    if (sqrtv < to_RR(0)) {
        return to_RR(-1);
    }

    // manually compute 2nd largest x-value (and corresponding y-value)
    xis[em-2] = to_RR(to_RR(sigma) * sqrt(sqrtv));
    rhos[em-2] = rho(xis[em-2], sigma);

    // compute the other x- and y-values iteratively from back to front
    for (int i = em-3; i > 0; i--) {
        sqrtv = m2 * log(area / to_RR(TruncToZZ(xis[i+1]) + to_ZZ(1)) + rho(xis[i+1], sigma));
        if (sqrtv < to_RR(0)) {
            return to_RR(-1);
        }
        xis[i] = to_RR(to_RR(sigma) * sqrt(sqrtv));
        rhos[i] = exp(-o2 * power(to_RR(xis[i] / sigma), 2));
    }

    // compute the y-value for x=0 and output it (as "quality" of the partition)
    rhos[0] = area / to_RR(TruncToZZ(xis[1]) + to_ZZ(1)) + rho(xis[1], sigma);
    return rhos[0];
}
