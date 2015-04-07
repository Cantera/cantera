/**
 * @file PrintCtrl.cpp
 *    Definitions for a simple class that augments the streams printing capabilities
 *   (see \ref Cantera::PrintCtrl).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/base/PrintCtrl.h"
#include "cantera/base/global.h"
#include <cmath>

using namespace std;

namespace Cantera
{

// Storage for the global crop flag
PrintCtrl::CROP_TYPE_GLOBAL PrintCtrl::GlobalCrop = GCT_NOPREF;

PrintCtrl::PrintCtrl(std::ostream& coutProxy, int Ndec,
                     CROP_TYPE ctlocal) :
    m_cout(coutProxy),
    m_Ndec(Ndec),
    m_precision(12),
    m_wMin(9),
    m_wMax(19),
    m_cropCntrl(ctlocal)
{
    warn_deprecated("class PrintCtrl");
}

void PrintCtrl::pr_de_c10(const double din, int p, const int wMin,
                          const int wMax)
{
    double d = cropAbs10(din, m_Ndec);
    pr_de(d, p, wMin, wMax);
}

void PrintCtrl::pr_de(const double d, int sigDigIn, const int wMinIn,
                      const int wMaxIn)
{
    int p = m_precision;
    if (sigDigIn != -1) {
        p = sigDigIn-1;
        if (p < 0) {
            p = 0;
        }
    }

    int wMin = m_wMin;
    if (wMinIn != -1) {
        wMin = wMinIn;
        if (wMin < 1) {
            wMin = 1;
        }
    }

    int wMax = m_wMax;
    if (wMaxIn != -1) {
        wMax = wMaxIn;
        if (wMax < 1) {
            wMax = 1;
        }
    }

    if (wMin > wMax) {
        wMax = wMin;
    }

    // Have to do the wMax ourselves, since C++ doesn't seem to
    // have a streams manipulator to do this !?!
    double dfabs = fabs(d);
    // This is the normal length assuming no sign and an 1.0E+04
    // formated exponented
    int requestedLength = 6 + p;
    if (d < 0.0) {
        requestedLength++;
    }
    if (dfabs < 9.9999999999E-99) {
        requestedLength++;
    }
    if (dfabs > 9.9999999999E99) {
        requestedLength++;
    }
    if (requestedLength > wMax) {
        p -= (requestedLength - wMax);
        if (p < 0) {
            p = 0;
        }
    }

    // Set to upper case and scientific notation
    m_cout.setf(ios_base::scientific | ios_base::uppercase);
    int wold = (int) m_cout.width(wMin);
    int pold = (int) m_cout.precision(p);

    m_cout << d;
    // Return the precision to the previous value;
    m_cout.precision(pold);
    m_cout.unsetf(ios_base::scientific);

    // Return width to original
    m_cout.width(wold);
}

double PrintCtrl::cropAbs10(const double d, int Ndec) const
{
    if (!doCrop()) {
        return d;
    }
    if (Ndec < -301 || Ndec > 301) {
        return d;
    }
    double dfabs = fabs(d);
    double  pdec = pow(10.0, (double) Ndec);
    if (dfabs < pdec) {
        return 0.0;
    }
    double dl10 = log10(dfabs);
    int N10 = (int) dl10;
    if (dl10 > -0.0) {
        N10 += 1;
    }
    int nsig = N10 - Ndec;
    return cropSigDigits(d, nsig);
}

double PrintCtrl::cropSigDigits(const double d, int nSig) const
{
    if (!doCrop()) {
        return d;
    }
    if (nSig <=0) {
        nSig = 1;
    }
    if (nSig >=9) {
        nSig = 9;
    }
    double sgn = 1.0;
    if (d < 0.0) {
        sgn = -1.0;
    }
    double dfabs = fabs(d);
    double dl10 = log10(dfabs);
    int N10 = (int) dl10;
    if (dl10 > -0.0) {
        N10 += 1;
    }
    int E10 = -N10 + nSig ;
    double pfabs = dfabs * pow(10.0, (double) E10);
    pfabs *= (1.0 + 1.0E-14);
    long int nfabs = (long int) pfabs;
    double remainder = pfabs - nfabs;
    if (remainder > 0.5) {
        nfabs++;
    }
    double paltabs = (double) nfabs;
    double daltabs = paltabs * pow(10.0, (double) -E10);
    return sgn * daltabs;
}

int PrintCtrl::setNdec(int Ndec)
{
    int nold = m_Ndec;
    m_Ndec = Ndec;
    return nold;
}

int PrintCtrl::setSigDigits(int nSigDigits)
{
    int nold = m_precision + 1;
    m_precision = nSigDigits - 1;
    if (m_precision < 0) {
        m_precision = 0;
    }
    return nold;
}

int PrintCtrl::setWmin(int wmin)
{
    int nold = m_wMin;
    m_wMin = wmin;
    return nold;
}

int PrintCtrl::setWmax(int wmax)
{
    int nold = m_wMax;
    m_wMax = wmax;
    return nold;
}

bool PrintCtrl::doCrop() const
{
    bool retn = ((m_cropCntrl == CT_ON) || (m_cropCntrl == CT_ON_GLOBALOBEY));
    if (m_cropCntrl ==  CT_ON_GLOBALOBEY) {
        if (GlobalCrop ==  GCT_NOCROP) {
            retn = false;
        }
    } else  if (m_cropCntrl == CT_OFF_GLOBALOBEY) {
        if (GlobalCrop == GCT_CROP) {
            retn = true;
        }
    }
    return retn;
}

void PrintCtrl:: setCropCntrl(CROP_TYPE ctlocal)
{
    m_cropCntrl = ctlocal;
}
}
