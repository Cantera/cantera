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
#include "cantera/numerics/eigen_dense.h"

namespace Cantera
{
/**
 * This class calculates the properties of electron in a gas.
 * @ingroup electron
 */
class Electron
{
public:
    Electron();

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
        return m_gridC[i];
    }

    //! Setup grid of electron energy.
    void setupGrid(size_t n, const double* eps);

    //! electron diffusivity
    virtual double electronDiffusivity() {
        throw NotImplementedError("Electron::electronDiffusivity");
    }

    //! electron mobility
    virtual double electronMobility() {
        throw NotImplementedError("Electron::electronMobility");
    }

    //! mean electron energy
    virtual double meanElectronEnergy() {
        throw NotImplementedError("Electron::meanElectronEnergy");
    }

    virtual double powerGain() {
        throw NotImplementedError("Electron::powerGain");
    }

    //! elastic power loss
    virtual double elasticPowerLoss() {
        throw NotImplementedError("Electron::elasticPowerLoss");
    }

    //! inelastic power loss
    virtual double inelasticPowerLoss() {
        throw NotImplementedError("Electron::inelasticPowerLoss");
    }

    //! total collision frequency
    virtual double totalCollisionFreq() {
        throw NotImplementedError("Electron::inelasticPowerLoss");
    }

    //! rate coefficient. [m^3/s]
    virtual double rateCoefficient(size_t k) {
        throw NotImplementedError("Electron::rateCoefficient");
    }

    //! inverse rate coefficient. [m^3/s]
    virtual double inverseRateCoefficient(size_t k) {
        throw NotImplementedError("Electron::inverseRateCoefficient");
    }

    //! net plasma production rates
    virtual void getNetPlasmaProductionRates(double* wdot) {
        throw NotImplementedError("Electron::getNetPlasmaProductionRates");
    }

    //! initialize Electron. Need to be called after adding all cross sections.
    void init(thermo_t* thermo);

    //! Reduced electric field
    double electricField() const {
        return m_E;
    }

    //! Set reduced electric field
    void setElectricField(double E) {
        if (m_E != E) {
            m_E = E;
            m_f0_ok = false;
        }
    }

    //! Electric field frequency
    double electricFieldFreq() const {
        return m_F;
    }

    //! Set electric field frequency
    void setElectricFieldFreq(double F) {
        if (m_F != F) {
            m_F = F;
            m_f0_ok = false;
        }
    }

    /**
     * Set the parameters for the Boltzmann solver
     * @param maxn Maximum number of iterations
     * @param rtol Relative tolerance. The iteration is stopped when the norm
     *        of the absolute difference between EEDFs is smaller than rtol.
     * @param delta0 Initial value of the iteration parameter. This parameter
     *        is adapted in succesive iterations to improve convergence.
     * @param m Reduction factor of error. The Richardson extrapolation attemps
     *        to reduce the error by a factor of m in each iteration. Larger m
     *        means faster convergence but also has higher risk of encountering
     *        numerical instabilities.
     * @param init_kTe Initial electron mean energy in [eV]. Assume
     *        initial EEDF to be Maxwell-Boltzmann distribution at init_kTe.
     * @param warn Flag of showing warning of insufficient cross section data.
     */
    void setBoltzmannSolver(size_t maxn, double rtol, double delta0,
                            double m, double init_kTe, bool warn) {
        m_maxn = maxn;
        m_rtol = rtol;
        m_delta0 = delta0;
        m_factorM = m;
        m_init_kTe = init_kTe;
        m_warn = warn;
    }

    //! Return electron temperature
    virtual double electronTemperature() {
        throw NotImplementedError("Electron::electronTemperature");
    }

    //! Set chemionization scattering-in rate
    //! Equal to the reaction rate divided by gas and electron number density
    virtual void setChemionScatRate(double rate) {
        throw NotImplementedError("Electron::setChemionScatRate");
    }

    //! Check that an array size is at least nSpecies()
    //! Throws an exception if kk is less than nSpecies(). Used before calls
    //! which take an array pointer.
    void checkSpeciesArraySize(size_t k) const;

    //! list of targets of electron collision
    std::vector<std::string> m_targets;

    //! list of kinds of electron collision
    std::vector<std::string> m_kinds;

    //! list of products of electron collision
    std::vector<std::string> m_products;

    //! list of mass ratio of electron to target species
    vector_fp m_massRatios;

    //! list of thresholds of electron collision
    vector_fp m_thresholds;

    //! cross sections
    std::vector<std::vector<std::vector<double>>> m_crossSections;

protected:
    //! Cached for saved calculations within each Electron.
    /*!
     *   For more information on how to use this, see examples within the source
     *   code and documentation for this within ValueCache class itself.
     */
    mutable ValueCache m_cache;

    //! Update temperature
    void update_T();

    //! Update composition
    void update_C();

    //! Calculate elastic cross section
    void calculateElasticCrossSection();

    //! number of cross section sets
    size_t m_ncs;

    //! Grid of electron energy (cell center) [eV]
    vector_fp m_gridC;

    //! Grid of electron energy (cell boundary i-1/2) [eV]
    vector_fp m_gridB;

    //! number of points for energy grid
    size_t m_points;

    //! Boltzmann constant times gas temperature
    double m_kT;

    //! reduced electric field
    double m_E;

    //! electric field freq
    double m_F;

    //! normalized electron energy distribution function
    Eigen::VectorXd m_f0;

    //! constant gamma
    double m_gamma;

    //! mole fractions of target
    vector_fp m_moleFractions;

    //! shift factor
    std::vector<int> m_shiftFactor;
    std::vector<int> m_inFactor;

    //! list elastic
    std::vector<size_t> m_kElastic;

    //! list inelastic
    std::vector<size_t> m_kInelastic;

    //! list effective
    std::vector<size_t> m_kEffective;

    //! list solo elastic
    std::vector<size_t> m_kSoloElastic;

    //! flag of electron energy distribution function
    bool m_f0_ok;

    //! pointer to the object representing the phase
    thermo_t* m_thermo;

    //! local gas composition
    compositionMap m_gasComposition;

    //! Maximum number of iterations
    size_t m_maxn;

    //! Relative tolerance
    double m_rtol;

    //! Initial value of the iteration parameter
    double m_delta0;

    //! Reduction factor of error
    double m_factorM;

    //! Initial electron mean energy
    double  m_init_kTe;

    //! Flag of warning of insufficient cross section data
    bool m_warn;

    //! Gas number density
    double m_N;
};

}

#endif
