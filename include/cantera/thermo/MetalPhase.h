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
    //! The electron cloud adds no volume of its own; any volume change of the
    //! electrode belongs to the host phase's equation of state.
    double molarVolume() const override {
        return 0.0;
    }
    double enthalpy_mole() const override {
        return 0.0;
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

    void getChemPotentials(span<double> mu) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu[n] = 0.0;
        }
    }

    void getEnthalpy_RT(span<double> hrt) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            hrt[n] = 0.0;
        }
    }

    void getEntropy_R(span<double> sr) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            sr[n] = 0.0;
        }
    }

    void getStandardChemPotentials(span<double> mu0) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu0[n] = 0.0;
        }
    }

    void getActivityConcentrations(span<double> c) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            c[n] = 1.0;
        }
    }

    void getPartialMolarEnthalpies(span<double> h) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            h[n] = 0.0;
        }
    }

    void getPartialMolarVolumes(span<double> vbar) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            vbar[n] = 0.0;
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
            warn_deprecated("MetalPhase::initThermo",
                "The 'density' field for 'electron-cloud' phases is deprecated "
                "and no longer affects the phase's properties; the molar volume "
                "is zero. To be removed after Cantera 4.0.");
            assignDensity(m_input.convert("density", "kg/m^3"));
        }
    }

private:
    double m_press;
};
}

#endif
