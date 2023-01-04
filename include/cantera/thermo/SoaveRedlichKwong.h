//! @file SoaveRedlichKwong.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOAVEREDLICHKWONG_H
#define CT_SOAVEREDLICHKWONG_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Soave-Redlich-Kwong equation of state
 *
 * @ingroup thermoprops
 */
class SoaveRedlichKwong : public MixtureFugacityTP
{
public:
    explicit SoaveRedlichKwong(const std::string& infile="",
                               const std::string& id="");

    virtual std::string type() const {
        return "Soave-Redlich-Kwong";
    }

public:
    int SoaveRedlichKwong::solveCubic(double T, double pres, double a, double b, double aAlpha,
                                  double Vroot[3]) const;

private:
    static const double omega_a;
    static const double omega_b;
    static const double omega_vc;

};
}

#endif
