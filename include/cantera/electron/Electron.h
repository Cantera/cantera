/**
 * @file Electron.h
 * Header file for class Electron.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRON_H
#define CT_ELECTRON_H

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/electron/ElectronCrossSection.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ValueCache.h"

namespace Cantera
{
/**
 * This class calculates the properties of electron in a gas.
 * @ingroup electron
 */
class Electron
{
public:
    Electron(); //!< Default constructor.

    virtual ~Electron();

    // Electron objects are not copyable or assignable
    Electron(const Electron&) = delete;
    Electron& operator=(const Electron&) = delete;

    //! Add a electron corss section to this Electron. Returns `true` if the electron cross section was
    //! successfully added, or `false` if the electron cross section was ignored.
    virtual bool addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs);

    //! number of cross section dataset
    size_t nElectronCrossSections() const {
        return m_ncs;
    }

    //! number of points of energy grid
    size_t nPoints() const {
        return m_points;
    }

    //! energy grid
    double grid(size_t i) const {
        return m_eps[i];
    }

    //! Setup grid of electron energy.
    void setupGrid(size_t n, const double* eps);

    //! Setup cross sections.
    void setupCrossSections();

    //! electron diffusivity
    //! @param N gas number density in SI
    double electronDiffusivity(double N);

    //! electron mobility
    //! @param N gas number density in SI
    double electronMobility(double N);

    void init(thermo_t* thermo);

    std::vector<std::string> m_electronCrossSectionTargets;
    std::vector<std::string> m_electronCrossSectionKinds;
    vector_fp m_massRatios;
    std::vector<std::vector<std::vector<double>>> m_electronCrossSectionData;

protected:
    //! Cached for saved calculations within each Electron.
    /*!
     *   For more information on how to use this, see examples within the source
     *   code and documentation for this within ValueCache class itself.
     */
    mutable ValueCache m_cache;

    //! calculate total cross section
    void calculateTotalCrossSection();

    //! calculate electron energy distribution function
    virtual void calculateDistributionFunction();

    //! normalized net production frequency
    //! Equation 10 of [1]
    //! vi / (N gamma)
    //! @param f0 normalized electron energy distribution function
    double netProductionFrequency(const vector_fp& f0);

    //! update temperature
    void update_T();

    //! update composition
    void update_C();

    //! number of cross section sets
    size_t m_ncs;

    //! grid of electron energy [eV]
    vector_fp m_eps;

    //! number of points for energy grid
    size_t m_points;

    //! Boltzmann constant times electron temperature
    double m_kTe;

    //! Boltzmann constant times gas temperature
    double m_kT;

    //! normalized electron energy distribution function
    vector_fp m_f0;

    //! the derivative of normalized electron energy distribution function
    vector_fp m_df0;

    //! constant gamma
    double m_gamma;

    //! mole fractions of target
    vector_fp m_moleFractions;

    //! flag of electron Cross Sections
    bool m_electronCrossSections_ok;

    //! flag of electron energy distribution function
    bool m_f0_ok;

    //! flag of total cross section
    bool m_totalCrossSection_ok;

    //! vector of electron cross section on the energy grid
    std::vector<vector_fp> m_electronCrossSections;

    //! vector of total electron cross section on the energy grid
    vector_fp m_totalCrossSection;

    //! vector of total attachment cross section on the energy grid
    vector_fp m_attachCrossSection;

    //! vector of total ionization cross section on the energy grid
    vector_fp m_ionizCrossSection;

    //! pointer to the object representing the phase
    thermo_t* m_thermo;

    //! local gas composition
    compositionMap m_gasComposition; 
};

}

#endif
