/**
 * @file PlasmaPhase.h
 * Header file for class PlasmaPhase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAPHASE_H
#define CT_PLASMAPHASE_H

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/plasma/ElectronCrossSection.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ValueCache.h"
#include "cantera/numerics/eigen_dense.h"

namespace Cantera
{
unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node);
/**
 * @defgroup electron
 * This class calculates the electron energy distribution function (EEDF) in a gas
 * by modeling collisions between electrons and other species represented by the
 * class ElectronCrossSection. EEDF is used to calculate reaction rate coefficient 
 * for plasma reaction and electron temperature for electron-temperature reaction
 * used in kinetics, and diffusivity/mobility of electron in transport.
 */

/*!
 * Class Plasma is the base class which manages the grid and grid cache,
 * cross-section data, and updating temperature and gas composition.
 * @ingroup electron
 */
class PlasmaPhase: public IdealGasPhase
{
public:
    PlasmaPhase();

    void addElectronCrossSections(const AnyValue& crossSections,
                                  const AnyValue& names);

    //! Add an electron cross section to this Electron. Returns `true` if the electron cross section was
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
        return m_gridCenter[i];
    }

    //! Setup grid of electron energy. n is the dimension of grid.
    //! eps is the vector of electron energy [eV].
    void setupGrid(size_t n, const double* eps);

    //! electron diffusivity
    virtual double electronDiffusivity() {
        throw NotImplementedError("PlasmaPhase::electronDiffusivity");
    }

    //! electron mobility
    virtual double electronMobility() {
        throw NotImplementedError("PlasmaPhase::electronMobility");
    }

    //! mean electron energy
    virtual double meanElectronEnergy() {
        throw NotImplementedError("PlasmaPhase::meanElectronEnergy");
    }

    virtual double powerGain() {
        throw NotImplementedError("PlasmaPhase::powerGain");
    }

    //! elastic power loss
    virtual double elasticPowerLoss() {
        throw NotImplementedError("PlasmaPhase::elasticPowerLoss");
    }

    //! inelastic power loss
    virtual double inelasticPowerLoss() {
        throw NotImplementedError("PlasmaPhase::inelasticPowerLoss");
    }

    //! total collision frequency
    virtual double totalCollisionFreq() {
        throw NotImplementedError("PlasmaPhase::totalCollisionFreq");
    }

    //! rate coefficient for the electron collision process. [m^3/s]
    virtual double rateCoefficient(size_t k) {
        throw NotImplementedError("PlasmaPhase::rateCoefficient");
    }

    //! reverse rate coefficient for the electron collision process. [m^3/s]
    virtual double reverseRateCoefficient(size_t k) {
        throw NotImplementedError("PlasmaPhase::reverseRateCoefficient");
    }

    //! initialize Plasma.
    virtual void initPlasma(const AnyMap& phaseNode, const AnyMap& rootNode);

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

    //! target of a specific process
    std::string target(size_t k) {
        return m_ecss[k]->target;
    }

    //! index of target
    size_t targetIndex(size_t k) {
        return m_kTargets[k];
    }

    //! kind of a specific process
    std::string kind(size_t k) {
        return m_ecss[k]->kind;
    }

    //! product of a specific process
    std::string product(size_t k) {
        return m_ecss[k]->product;
    }

    //! threshold of a specific process
    double threshold(size_t k) {
        return m_ecss[k]->threshold;
    }

    //! scattering-in factor
    int inFactor(size_t k) {
        return m_inFactor[k];
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
     */
    void setupBoltzmannSolver(size_t maxn, double rtol,
                              double delta0, double m) {
        m_maxn = maxn;
        m_rtol = rtol;
        m_delta0 = delta0;
        m_factorM = m;
    }

    /**
     * Set the threshold of mole fraction for showing warning of
     * insufficient cross-section data. The warning will show if
     * any substantial species whose mole fraction is higher than
     * the threshold lack the cross-section data.
     */
    void setMoleFractionThreshold(double fraction) {
        m_moleFractionThreshold = fraction;
    }

    /**
     * Set initial mean electron energy in [eV]. Assume initial
     * EEDF to be Maxwell-Boltzmann distribution at init_kTe.
     */
    void setInitialMeanElectronEnergy(double init_kTe) {
        m_init_kTe = init_kTe;
    }

    //! Set chemionization scattering-in rate
    //! Equal to the reaction rate divided by gas and electron number density
    virtual void setChemionScatRate(double rate) {
        throw NotImplementedError("PlasmaPhase::setChemionScatRate");
    }

    //! Overload to signal updating electron energy density function.
    virtual void setTemperature(const double temp);

protected:
    // set grid cache
    void setGridCache();

    //! Signal updating electron energy density function
    virtual void compositionChanged();

    //! check gas comsition for any substantial species without
    //! the cross-section data.
    void checkSpeciesNoCrossSection();

    //! Calculate elastic cross section
    void calculateElasticCrossSection();

    //! number of cross section sets
    size_t m_ncs;

    //! array of cross-section object
    std::vector<shared_ptr<ElectronCrossSection>> m_ecss;

    //! Grid of electron energy (cell center) [eV]
    vector_fp m_gridCenter;

    //! Grid of electron energy (cell boundary i-1/2) [eV]
    vector_fp m_gridEdge;

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

    //! Mole fractions of target species of each collision process
    vector_fp m_moleFractions;

    //! shift factor
    std::vector<int> m_shiftFactor;
    std::vector<int> m_inFactor;

    //! Indices of elastic collisions in m_crossSections
    std::vector<size_t> m_kElastic;

    //! Indices of inelastic collisions in m_crossSections
    std::vector<size_t> m_kInelastic;

    //! Indices of effective collisions in m_crossSections
    std::vector<size_t> m_kEffective;

    //! flag of electron energy distribution function
    bool m_f0_ok;

    //! Maximum number of iterations
    size_t m_maxn;

    //! Relative tolerance of electron energy distribution function
    //! for solving Boltzmann equation
    double m_rtol;

    //! Initial value of the iteration parameter
    double m_delta0;

    //! Reduction factor of error
    double m_factorM;

    //! Initial electron mean energy
    double  m_init_kTe;

    //! The threshold of mole fraction for showing warning of
    //! insufficient cross-section data.
    double m_moleFractionThreshold;

    //! Gas number density
    double m_N;

    //! Location of cell j for grid cache
    std::vector<std::vector<size_t>> m_j;

    //! Location of cell i for grid cache
    std::vector<std::vector<size_t>> m_i;

    //! Cross section at the boundaries of the overlap of cell i and j
    std::vector<std::vector<vector_fp>> m_sigma;

    //! The energy boundaries of the overlap of cell i and j
    std::vector<std::vector<vector_fp>> m_eps;

    //! cross section data. m_crossSections[i][j][k] where i is the specific process,
    //! j=0 is the vector of electron energy [eV], j=1 is the vector of cross section, and
    // k is the index of vector.
    std::vector<std::vector<std::vector<double>>> m_crossSections;

    //! list of target species indices
    std::vector<size_t> m_kTargets;

    //! Indices of species which has no cross-section data
    std::vector<size_t> m_kOthers;

    //! list of product index
    std::vector<std::vector<size_t>> m_kProducts;
};

}

#endif
