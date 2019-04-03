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

    size_t nElectronCrossSections() const {
        return m_ncs;
    }

    size_t nPoints() const {
        return m_points;
    }

    double grid(size_t i) const {
        return m_eps[i];
    }

    // Setup grid of electron energy.
    void setupGrid(size_t n, const double* eps);

    // Setup cross sections.
    void setupCrossSections();

    double electronDiffusivity(double N);
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

    void calculateTotalCrossSection();

    virtual void calculateDistributionFunction();

    // update temperature
    void update_T();

    // update composition
    void update_C();

    // Number of cross section sets
    size_t m_ncs;

    // Grid of electron energy [eV]
    vector_fp m_eps;

    // Number of points for energy grid
    size_t m_points;

    // Boltzmann constant times electron temperature
    double m_kTe;

    // Boltzmann constant times gas temperature
    double m_kT;

    // Electron energy distribution function
    vector_fp m_f0;
    vector_fp m_df0;

    double m_gamma;

    vector_fp m_moleFractions;
    double m_N;

    // Flag
    bool m_electronCrossSections_ok;
    bool m_f0_ok;

    // Vector of electron cross section on the energy grid
    std::vector<vector_fp> m_electronCrossSections;

    // Vector of total electron cross section on the energy grid
    vector_fp m_totalCrossSection;

    //! pointer to the object representing the phase
    thermo_t* m_thermo;

    //! stroe a local gas composition
    compositionMap m_gasComposition; 
};

//! typedef for the Electron class
typedef Electron electron_t;

}

#endif
