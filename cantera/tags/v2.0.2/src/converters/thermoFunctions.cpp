/**
 *  @file thermoFunctions.cpp
 *  File containing thermo evaluation functions for NASA polynomials,
 *  which are used in testing the interpolations.
 */

// Copyright 2001  California Institute of Technology

#include <math.h>
#include "thermoFunctions.h"
#include <iostream>
using namespace std;

namespace ckr
{

/**
 *  non-dimensional heat capacity (\f$ C_p/R \f$) at constant P for
 *  one species @param t temperature @param s species object
 */
double cp(double t, const Species& s)
{
    if (s.thermoFormatType == 1) {
        const vector_fp* cpc;
        int ireg = -1;
        for (int i = 0; i < s.nTempRegions; i++) {
            if (t <= s.maxTemps[i]) {
                ireg = i;
                break;
            }
        }
        cpc = s.region_coeffs[ireg];
        const vector_fp& c = *cpc;
        double cp0r = c[0]/(t*t) + c[1]/t + c[2] + c[3]*t + c[4]*t*t
                      + c[5]*t*t*t + c[6]*t*t*t*t;
        return cp0r;
    }
    const vector_fp* cpc;
    if (t > s.tmid) {
        cpc = &s.highCoeffs;
    } else {
        cpc = &s.lowCoeffs;
    }
    const vector_fp& c = *cpc;
    double cp0r = c[0] + c[1]*t + c[2]*t*t + c[3]*t*t*t + c[4]*t*t*t*t;
    return cp0r;
}



/**
 *  enthalpy in Kelvin (\f$ H/R \f$) for
 *  one species. @param t temperature @param s species object
 */
double enthalpy(double t, const Species& s)
{
    if (s.thermoFormatType == 1) {
        const vector_fp* cpc;
        int ireg = -1;
        for (int i = 0; i < s.nTempRegions; i++) {
            if (t <= s.maxTemps[i]) {
                ireg = i;
                break;
            }
        }
        cpc = s.region_coeffs[ireg];
        const vector_fp& c = *cpc;
        double h0rt = -c[0]/(t*t) + c[1]*log(t)/t
                      + c[2] + 0.5*c[3]*t + c[4]*t*t/3.0 + 0.25*c[5]*t*t*t
                      + 0.2*c[6]*t*t*t*t + c[7]/t;
        return t*h0rt;
    }
    const vector_fp* cp;
    if (t > s.tmid) {
        cp = &s.highCoeffs;
    } else {
        cp = &s.lowCoeffs;
    }
    const vector_fp& c = *cp;
    double h0rt = c[0] + 0.5*c[1]*t + c[2]*t*t/3.0 + 0.25*c[3]*t*t*t
                  + 0.2*c[4]*t*t*t*t + c[5]/t;
    return t*h0rt;
}


/**
 *  non-dimensional entropy (\f$ S/R \f$) for
 *  one species @param t temperature @param s species object
 */
double entropy(double t, const Species& s)
{
    if (s.thermoFormatType == 1) {
        const vector_fp* cpc;
        int ireg = -1;
        for (int i = 0; i < s.nTempRegions; i++) {
            if (t <= s.maxTemps[i]) {
                ireg = i;
                break;
            }
        }
        cpc = s.region_coeffs[ireg];
        const vector_fp& c = *cpc;
        double s0r =  -0.5*c[0]/(t*t) - c[1]/t
                      + c[2]*log(t) + c[3]*t + 0.5*c[4]*t*t + c[5]*t*t*t/3.0
                      + 0.25*c[6]*t*t*t*t + c[8];
        return t*s0r;
    }
    const vector_fp* cp;
    if (t > s.tmid) {
        cp = &s.highCoeffs;
    } else {
        cp = &s.lowCoeffs;
    }
    const vector_fp& c = *cp;
    double s0r = c[0]*log(t) + c[1]*t + 0.5*c[2]*t*t + c[3]*t*t*t/3.0
                 + 0.25*c[4]*t*t*t*t + c[6];
    return t*s0r;
}

/**
 *  Gibbs function in Kelvin (\f$ G/R \f$) for
 *  one species. @param t temperature @param s species object
 */
double gibbs(double t, const Species& s)
{
    if (s.thermoFormatType == 1) {
        double s0r = entropy(t, s);
        double h0r = enthalpy(t, s);
        return (h0r - s0r * t);
    }
    const vector_fp* cp;
    if (t > s.tmid) {
        cp = &s.highCoeffs;
    } else {
        cp = &s.lowCoeffs;
    }
    const vector_fp& c = *cp;
    double h0rt = c[0] + 0.5*c[1]*t + c[2]*t*t/3.0 + 0.25*c[3]*t*t*t
                  + 0.2*c[4]*t*t*t*t + c[5]/t;
    double s0r = c[0]*log(t) + c[1]*t + 0.5*c[2]*t*t + c[3]*t*t*t/3.0
                 + 0.25*c[4]*t*t*t*t + c[6];
    return t*(h0rt - s0r);
}

}

