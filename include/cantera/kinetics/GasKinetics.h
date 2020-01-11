/**
 * @file GasKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include "BulkKinetics.h"
#include "ThirdBodyCalc.h"
#include "FalloffMgr.h"
#include "Reaction.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary gas-phase chemistry. This kinetics manager
 * implements standard mass-action reaction rate expressions for low-density
 * gases.
 * @ingroup kinetics
 */
class GasKinetics : public BulkKinetics
{
public:
    //! @name Constructors and General Information
    //! @{

    //! Constructor.
    /*!
     *  @param thermo  Pointer to the gas ThermoPhase (optional)
     */
    GasKinetics(thermo_t* thermo = 0);

    virtual std::string kineticsType() const {
        return "Gas";
    }

    //! @}
    //! @name Reaction Rates Of Progress
    //! @{

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(doublereal* kfwd);

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{
    virtual void init();
    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
    virtual void invalidateCache();
    //@}

    void updateROP();

    //! Update temperature-dependent portions of reaction rates and falloff
    //! functions.
    virtual void update_rates_T();

    //! Update properties that depend on concentrations.
    //! Currently the enhanced collision partner concentrations are updated
    //! here, as well as the pressure-dependent portion of P-log and Chebyshev
    //! reactions.
    virtual void update_rates_C();

    //! enable kinetics to use class PlasmaElectron
    virtual void enableElectron(bool enable) {
        m_do_electron = enable;
    }

    //! Set the value of electron temperature in Kelvin. The electron temperature is equal
    //! to gas temperature by default, and can be set manually to a different temperature.
    //! However, if m_do_electron is set to true, class PlasmaElectron is used to calculate
    //! the value of electron temperature.
    virtual void setElectronTemperature(double Te) {
        m_Te_fix = Te;
    }

protected:
    //! Reaction index of each falloff reaction
    std::vector<size_t> m_fallindx;

    //! Map of reaction index to falloff reaction index (i.e indices in
    //! #m_falloff_low_rates and #m_falloff_high_rates)
    std::map<size_t, size_t> m_rfallindx;

    //! Rate expressions for falloff reactions at the low-pressure limit
    Rate1<Arrhenius> m_falloff_low_rates;

    //! Rate expressions for falloff reactions at the high-pressure limit
    Rate1<Arrhenius> m_falloff_high_rates;

    //! Rate expressions for electron reactions
    Rate1<ElectronArrhenius> m_electron_temperature_rates;

    //! Rate expressions for plasma reactions
    Rate1<PlasmaRate> m_plasma_rates;

    FalloffMgr m_falloffn;

    ThirdBodyCalc m_3b_concm;
    ThirdBodyCalc m_falloff_concm;

    Rate1<Plog> m_plog_rates;
    Rate1<ChebyshevRate> m_cheb_rates;

    //! @name Reaction rate data
    //!@{
    doublereal m_logp_ref;
    doublereal m_logc_ref;
    doublereal m_logStandConc;
    vector_fp m_rfn_low;
    vector_fp m_rfn_high;

    doublereal m_pres; //!< Last pressure at which rates were evaluated
    vector_fp falloff_work;
    vector_fp concm_3b_values;
    vector_fp concm_falloff_values;
    //!@}

    //! Electron temperature for electron-temperature reactions
    double m_temp_e;

    //! flag of enabling electrons
    bool m_do_electron;

    //! fixed value of electron temperature
    double m_Te_fix;

    void processFalloffReactions();
    void processPlasmaReactions();

    void addElectronTemperatureReaction(ElectronTemperatureReaction& r);
    void addPlasmaReaction(PlasmaReaction& r);
    void addThreeBodyReaction(ThreeBodyReaction& r);
    void addFalloffReaction(FalloffReaction& r);
    void addPlogReaction(PlogReaction& r);
    void addChebyshevReaction(ChebyshevReaction& r);

    void modifyElectronTemperatureReaction(size_t i, ElectronTemperatureReaction& r);
    void modifyPlasmaReaction(size_t i, PlasmaReaction& r);
    void modifyThreeBodyReaction(size_t i, ThreeBodyReaction& r);
    void modifyFalloffReaction(size_t i, FalloffReaction& r);
    void modifyPlogReaction(size_t i, PlogReaction& r);
    void modifyChebyshevReaction(size_t i, ChebyshevReaction& r);

    //! Update the equilibrium constants in molar units.
    void updateKc();
};

}

#endif
