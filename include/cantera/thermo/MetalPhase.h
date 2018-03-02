/**
 *  @file MetalPhase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_METALPHASE_H
#define CT_METALPHASE_H

#include "ThermoPhase.h"
#include "cantera/base/ctml.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 *
 * Class MetalPhase represents electrons in a metal.
 */
class MetalPhase : public ThermoPhase
{
public:
    MetalPhase() {}

    // Overloaded methods of class ThermoPhase

    virtual std::string type() const {
        return "Metal";
    }

    virtual double enthalpy_mole() const {
        return 0.0;
    }
    virtual double intEnergy_mole() const {
        return 0.0;
    }
    virtual double entropy_mole() const {
        return 0.0;
    }
    virtual double gibbs_mole() const {
        return 0.0;
    }
    virtual double cp_mole() const {
        return 0.0;
    }
    virtual double cv_mole() const {
        return 0.0;
    }

    virtual void setPressure(double pres) {
        m_press = pres;
    }
    virtual double pressure() const {
        return m_press;
    }

    virtual void getChemPotentials(double* mu) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu[n] = 0.0;
        }
    }

    virtual void getEnthalpy_RT(double* hrt) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            hrt[n] = 0.0;
        }
    }

    virtual void getEntropy_R(double* sr) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            sr[n] = 0.0;
        }
    }

    virtual void getStandardChemPotentials(double* mu0) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu0[n] = 0.0;
        }
    }

    virtual void getActivityConcentrations(double* c) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            c[n] = 1.0;
        }
    }

    virtual double standardConcentration(size_t k=0) const {
        return 1.0;
    }

    virtual double logStandardConc(size_t k=0) const {
        return 0.0;
    }

    virtual void setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","Metal");
        double rho = getFloat(eosdata, "density", "density");
        setDensity(rho);
    }

private:
    double m_press;
};
}

#endif
