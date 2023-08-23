/**
 *  @file MetalPhase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_METALPHASE_H
#define CT_METALPHASE_H

#include "ThermoPhase.h"

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

    string type() const override {
        return "electron-cloud";
    }

    bool isCompressible() const override {
        return false;
    }

    double enthalpy_mole() const override {
        return 0.0;
    }
    double intEnergy_mole() const override {
        return - pressure() * molarVolume();
    }
    double entropy_mole() const override {
        return 0.0;
    }
    double gibbs_mole() const override {
        return 0.0;
    }
    double cp_mole() const override {
        return 0.0;
    }
    double cv_mole() const override {
        return 0.0;
    }

    void setPressure(double pres) override {
        m_press = pres;
    }
    double pressure() const override {
        return m_press;
    }

    void getChemPotentials(double* mu) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu[n] = 0.0;
        }
    }

    void getEnthalpy_RT(double* hrt) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            hrt[n] = 0.0;
        }
    }

    void getEntropy_R(double* sr) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            sr[n] = 0.0;
        }
    }

    void getStandardChemPotentials(double* mu0) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu0[n] = 0.0;
        }
    }

    void getActivityConcentrations(double* c) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            c[n] = 1.0;
        }
    }
    void getPartialMolarEnthalpies(double *h) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            h[n] = 0.0;
        }
    }

    Units standardConcentrationUnits() const override {
        return Units(1.0);
    }

    double standardConcentration(size_t k=0) const override {
        return 1.0;
    }

    double logStandardConc(size_t k=0) const override {
        return 0.0;
    }

    void initThermo() override {
        if (m_input.hasKey("density")) {
            assignDensity(m_input.convert("density", "kg/m^3"));
        }
    }

    void getParameters(AnyMap& phaseNode) const override {
        ThermoPhase::getParameters(phaseNode);
        phaseNode["density"].setQuantity(density(), "kg/m^3");
    }

private:
    double m_press;
};
}

#endif
