/**
 *  @file MetalPhase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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

    virtual bool isCompressible() const {
        return false;
    }

    virtual doublereal enthalpy_mole() const {
        return 0.0;
    }
    virtual doublereal intEnergy_mole() const {
        return 0.0;
    }
    virtual doublereal entropy_mole() const {
        return 0.0;
    }
    virtual doublereal gibbs_mole() const {
        return 0.0;
    }
    virtual doublereal cp_mole() const {
        return 0.0;
    }
    virtual doublereal cv_mole() const {
        return 0.0;
    }

    virtual void setPressure(doublereal pres) {
        m_press = pres;
    }
    virtual doublereal pressure() const {
        return m_press;
    }

    virtual void getChemPotentials(doublereal* mu) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu[n] = 0.0;
        }
    }

    virtual void getEnthalpy_RT(doublereal* hrt) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            hrt[n] = 0.0;
        }
    }

    virtual void getEntropy_R(doublereal* sr) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            sr[n] = 0.0;
        }
    }

    virtual void getStandardChemPotentials(doublereal* mu0) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu0[n] = 0.0;
        }
    }

    virtual void getActivityConcentrations(doublereal* c) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            c[n] = 1.0;
        }
    }
    virtual void getPartialMolarEnthalpies(doublereal *h) const {
        for (size_t n = 0; n < nSpecies(); n++) {
            h[n] = 0.0;
        }
    }

    virtual Units standardConcentrationUnits() const {
        return Units(1.0);
    }

    virtual doublereal standardConcentration(size_t k=0) const {
        return 1.0;
    }

    virtual doublereal logStandardConc(size_t k=0) const {
        return 0.0;
    }

    virtual void initThermo() {
        if (m_input.hasKey("density")) {
            assignDensity(m_input.convert("density", "kg/m^3"));
        }
    }

    virtual void setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","Metal");
        doublereal rho = getFloat(eosdata, "density", "density");
        assignDensity(rho);
    }

private:
    doublereal m_press;
};
}

#endif
