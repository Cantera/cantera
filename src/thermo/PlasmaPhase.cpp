//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include <boost/math/special_functions/gamma.hpp>
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"
#include "cantera/numerics/funcs.h"

namespace Cantera {

PlasmaPhase::PlasmaPhase(const std::string& inputFile, const std::string& id_)
    : m_isotropicShapeFactor(2.0)
    , m_nPoints(1001)
    , m_electronSpeciesIndex(npos)
    , m_distributionType("isotropic")
{
    initThermoFile(inputFile, id_);

    // initial grid
    m_electronEnergyLevels = Eigen::ArrayXd::LinSpaced(m_nPoints, 0.0, 1.0);

    // initial electron temperature
    setElectronTemperature(temperature());
}

void PlasmaPhase::updateElectronEnergyDistribution()
{
    if (m_distributionType == "user-specified") {
        throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
            "Invalid for user-specified electron energy distribution.");
    } else if (m_distributionType == "isotropic") {
        setIsotropicElectronEnergyDistribution();
    }
}

void PlasmaPhase::setElectronEnergyDistributionType(const std::string& type)
{
    if (type == "user-specified" ||
        type == "isotropic") {
        m_distributionType = type;
    } else {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistributionType",
            "Unknown type for electron energy distribution.");
    }
}

void PlasmaPhase::setIsotropicElectronEnergyDistribution()
{
    m_electronEnergyDist.resize(m_nPoints);
    double x = m_isotropicShapeFactor;
    double gamma1 = boost::math::tgamma(3.0 / 2.0 * x);
    double gamma2 = boost::math::tgamma(5.0 / 2.0 * x);
    double c1 = x * std::pow(gamma2, 1.5) / std::pow(gamma1, 2.5);
    double c2 = x * std::pow(gamma2 / gamma1, x);
    m_electronEnergyDist =
        c1 * m_electronEnergyLevels.sqrt() /
        std::pow(m_meanElectronEnergy, 1.5) *
        (-c2 * (m_electronEnergyLevels /
        m_meanElectronEnergy).pow(x)).exp();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    m_electronTemp = Te;
    m_meanElectronEnergy = 3.0 / 2.0 * electronTemperature() *
                           Boltzmann / ElectronCharge;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setElectronEnergyLevels(const vector_fp& levels)
{
    m_nPoints = levels.size();
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::ArrayXd>(levels.data(), levels.size());
    updateElectronEnergyDistribution();
}

void PlasmaPhase::getElectronEnergyLevels(vector_fp& levels) const
{
    levels.resize(m_nPoints);
    Eigen::Map<Eigen::ArrayXd>(levels.data(), m_nPoints) = m_electronEnergyLevels;
}

void PlasmaPhase::setElectronEnergyDistribution(const vector_fp& levels,
                                                const vector_fp& distrb)
{
    m_distributionType = "user-specified";
    if (levels.size() != distrb.size()) {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistribution",
                           "Vector lengths need to be the same.");
    }
    m_nPoints = levels.size();
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::ArrayXd>(levels.data(), levels.size());
    m_electronEnergyDist =
        Eigen::Map<const Eigen::VectorXd>(distrb.data(), distrb.size());
    // calculate mean electron energy and electron temperature
    Eigen::ArrayXd eps52 = m_electronEnergyLevels.pow(5./2.);
    m_meanElectronEnergy =
        2.0 / 5.0 * trapezoidal(m_electronEnergyDist, eps52);
    Phase::setElectronTemperature(
        2.0 / 3.0 * m_meanElectronEnergy * Avogadro *
        ElectronCharge / GasConstant);
}

void PlasmaPhase::getElectronEnergyDistribution(vector_fp& distrb) const
{
    distrb.resize(m_nPoints);
    Eigen::Map<Eigen::VectorXd>(distrb.data(), m_nPoints) = m_electronEnergyDistrb;
}

void PlasmaPhase::setIsotropicShapeFactor(double x) {
    m_isotropicShapeFactor = x;
    setIsotropicElectronEnergyDistribution();
}

void PlasmaPhase::getParameters(AnyMap& phaseNode) const
{
    IdealGasPhase::getParameters(phaseNode);
    AnyMap eedf;
    eedf["type"] = m_distributionType;
    eedf["shape-factor"] = AnyValue(m_isotropicShapeFactor);
    phaseNode["electron-energy-distribution"] = std::move(eedf);
}

void PlasmaPhase::setParameters(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    IdealGasPhase::setParameters(phaseNode, rootNode);
    if (phaseNode.hasKey("electron-energy-distribution")) {
        const AnyMap eedf = phaseNode["electron-energy-distribution"].as<AnyMap>();
        m_distributionType = eedf["type"].asString();
        m_isotropicShapeFactor = eedf["shape-factor"].asDouble();
    }
}

bool PlasmaPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = IdealGasPhase::addSpecies(spec);
    size_t k = m_kk - 1;

    if (spec->composition.find("E") != spec->composition.end() &&
        spec->composition.size() == 1 &&
        spec->composition["E"] == 1) {
        if (m_electronSpeciesIndex == npos) {
            m_electronSpeciesIndex = k;
        } else {
            throw CanteraError("PlasmaPhase::addSpecies",
                               "Cannot add species, {}. "
                               "Only one electron species is allowed.", spec->name);
        }
    }
    return added;
}

void PlasmaPhase::initThermo()
{
    IdealGasPhase::initThermo();
    // check electron species
    if (m_electronSpeciesIndex == npos) {
        throw CanteraError("PlasmaPhase::initThermo",
                           "No electron species found.");
    }
}

void PlasmaPhase::updateThermo() const
{
    IdealGasPhase::updateThermo();
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tempNow = temperature();
    double electronTempNow = electronTemperature();
    size_t k = m_electronSpeciesIndex;
    // If the electron temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tempNow || cached.state2 != electronTempNow) {
        m_spthermo.update_single(k, electronTemperature(),
                &m_cp0_R[k], &m_h0_RT[k], &m_s0_R[k]);
        cached.state1 = tempNow;
        cached.state2 = electronTempNow;
    }
    // update the species Gibbs functions
    m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
}

}
