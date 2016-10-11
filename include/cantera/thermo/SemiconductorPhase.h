/**
 *  @file SemiconductorPhase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_SEMICONDPHASE_H
#define CT_SEMICONDPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"
#include "cantera/base/ctml.h"

namespace Cantera
{

const int cElectron = 0;
const int cHole = 1;

/**
 * @ingroup thermoprops
 *
 * Class SemiconductorPhase represents electrons and holes in a semiconductor.
 * @deprecated Broken and unused. To be removed after Cantera 2.3.
 *
 */
class SemiconductorPhase : public ThermoPhase
{
public:
    SemiconductorPhase() {
        warn_deprecated("class SemiconductorPhase",
                        "To be removed after Cantera 2.3.");
    }
    SemiconductorPhase(std::string infile, std::string id="");

    SemiconductorPhase(const SemiconductorPhase& right) {
        *this = right;
    }

    SemiconductorPhase& operator=(const SemiconductorPhase& right) {
        if (&right != this) {
            ThermoPhase::operator=(right);
            m_press = right.m_press;
        }
        return *this;
    }

    virtual ThermoPhase* duplMyselfAsThermoPhase() const {
        SemiconductorPhase* idg = new SemiconductorPhase(*this);
        return (ThermoPhase*) idg;
    }

    // Overloaded methods of class ThermoPhase

    virtual int eosType() const {
        warn_deprecated("SemiconductorPhase::eosType",
                        "To be removed after Cantera 2.3.");
        return cSemiconductor;
    }
    virtual std::string type() const {
        return "Semiconductor";
    }

    virtual void setPressure(doublereal pres) {
        m_press = pres;
    }
    virtual doublereal pressure() const {
        return m_press;
    }

    virtual void setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","Semiconductor");
        doublereal rho = getFloat(eosdata, "density", "-");
        setDensity(rho);
        m_bandgap = getFloat(eosdata, "bandgap", "-");
        doublereal e_mass = getFloat(eosdata, "electron_mass", "-");
        doublereal h_mass = getFloat(eosdata, "hole_mass", "-");
        doublereal e_donor = getFloat(eosdata, "donor_energy", "-");
        doublereal n_donor = getFloat(eosdata, "donor_concentration", "-");
        doublereal e_acceptor = getFloat(eosdata, "acceptor_energy", "-");
        doublereal n_acceptor = getFloat(eosdata, "acceptor_concentration", "-");
        setEffectiveMasses(e_mass, h_mass);
        setDonorDoping(n_donor, e_donor);
        setAcceptorDoping(n_acceptor, e_acceptor);
    }

    void setEffectiveMasses(doublereal e_mass, doublereal h_mass) {
        m_emass = e_mass;
        m_hmass = h_mass;
    }

    void setDonorDoping(doublereal n_donor, doublereal e_donor) {
        m_ndonor = n_donor;
        m_edonor = e_donor;
    }

    void setAcceptorDoping(doublereal n_acceptor, doublereal e_acceptor) {
        m_nacceptor = n_acceptor;
        m_eacceptor = e_acceptor;
    }

    doublereal effectiveMass_e() const {
        return m_emass;
    }

    doublereal effectiveMass_h() const {
        return m_hmass;
    }

    doublereal fermiLevel() const {
        return m_fermi_level;
    }

    virtual void getChemPotentials(doublereal* mu) const;
    doublereal nc() const;
    doublereal nv() const;

    /*!
     * Energy at the top of the conduction band. By default, energies are
     * referenced to this energy, and so this function simply returns zero.
     */
    doublereal ec() const;
    doublereal ev() const;
    doublereal bandgap() const {
        return m_bandgap;
    }

private:
    doublereal m_press;
    doublereal m_emass;
    doublereal m_hmass;
    doublereal m_ndonor;
    doublereal m_edonor;
    doublereal m_nacceptor;
    doublereal m_eacceptor;
    doublereal m_fermi_level;
    doublereal m_bandgap;
    mutable vector_fp m_work;

    void initLengths();
};
}

#endif
