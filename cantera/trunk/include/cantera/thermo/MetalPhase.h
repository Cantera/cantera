/**
 *  @file MetalPhase.h
 */

//  Copyright 2003 California Institute of Technology

#ifndef CT_METALPHASE_H
#define CT_METALPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 *
 * Class MetalPhase represents electrons in a metal.
 *
 */
class MetalPhase : public ThermoPhase
{

public:

    MetalPhase() {}

    MetalPhase(const MetalPhase& right) {
        *this = operator=(right);
    }

    MetalPhase& operator=(const MetalPhase& right) {
        if (&right != this) {
            ThermoPhase::operator=(right);
            m_press = right.m_press;
        }
        return *this;
    }

    //! Duplicator
    virtual ThermoPhase* duplMyselfAsThermoPhase() const {
        MetalPhase* idg = new MetalPhase(*this);
        return (ThermoPhase*) idg;
    }

    // Overloaded methods of class ThermoPhase

    virtual int eosType() const {
        return cMetal;
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
    virtual doublereal  pressure() const {
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

    virtual doublereal standardConcentration(size_t k=0) const {
        return 1.0;
    }

    virtual doublereal logStandardConc(size_t k=0) const {
        return 0.0;
    }

    virtual void setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","Metal");
        doublereal rho = ctml::getFloat(eosdata, "density", "density");
        setDensity(rho);
    }

protected:

private:
    doublereal m_press;
};
}

#endif





