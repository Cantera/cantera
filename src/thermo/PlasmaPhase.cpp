//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/polymorphic_pointer_cast.hpp>
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"
#include "cantera/numerics/funcs.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"

namespace Cantera {

namespace {
    const double gamma = sqrt(2 * ElectronCharge / ElectronMass);
}

PlasmaPhase::PlasmaPhase(const string& inputFile, const string& id_)
{
    initThermoFile(inputFile, id_);

    // initial grid
    m_electronEnergyLevels = Eigen::ArrayXd::LinSpaced(m_nPoints, 0.0, 1.0);

    // initial electron temperature
    setElectronTemperature(temperature());
}

void PlasmaPhase::updateElectronEnergyDistribution()
{
    if (m_distributionType == "discretized") {
        throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
            "Invalid for discretized electron energy distribution.");
    } else if (m_distributionType == "isotropic") {
        setIsotropicElectronEnergyDistribution();
    }
    updateElectronEnergyDistDifference();
    electronEnergyDistributionChanged();
}

void PlasmaPhase::normalizeElectronEnergyDistribution() {
    Eigen::ArrayXd eps32 = m_electronEnergyLevels.pow(3./2.);
    double norm = 2./3. * numericalQuadrature(m_quadratureMethod,
                                              m_electronEnergyDist, eps32);
    if (norm < 0.0) {
        throw CanteraError("PlasmaPhase::normalizeElectronEnergyDistribution",
                           "The norm is negative. This might be caused by bad "
                           "electron energy distribution");
    }
    m_electronEnergyDist /= norm;
}

void PlasmaPhase::setElectronEnergyDistributionType(const string& type)
{
    if (type == "discretized" ||
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
        std::pow(meanElectronEnergy(), 1.5) *
        (-c2 * (m_electronEnergyLevels /
        meanElectronEnergy()).pow(x)).exp();
    checkElectronEnergyDistribution();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    m_electronTemp = Te;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setMeanElectronEnergy(double energy) {
    m_electronTemp = 2.0 / 3.0 * energy * ElectronCharge / Boltzmann;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setElectronEnergyLevels(const double* levels, size_t length)
{
    m_nPoints = length;
    m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(levels, length);
    checkElectronEnergyLevels();
    electronEnergyLevelChanged();
    updateElectronEnergyDistribution();
    m_interp_cs.resize(m_nPoints);
    // The cross sections are interpolated on the energy levels
    if (nCollisions() > 0) {
        updateInterpolatedCrossSections();
    }
}

void PlasmaPhase::electronEnergyDistributionChanged()
{
    m_distNum++;
}

void PlasmaPhase::electronEnergyLevelChanged()
{
    m_levelNum++;
}

void PlasmaPhase::interpolatedCrossSectionsChanged()
{
    m_csNum++;
}

void PlasmaPhase::checkElectronEnergyLevels() const
{
    Eigen::ArrayXd h = m_electronEnergyLevels.tail(m_nPoints - 1) -
                       m_electronEnergyLevels.head(m_nPoints - 1);
    if (m_electronEnergyLevels[0] < 0.0 || (h <= 0.0).any()) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
            "Values of electron energy levels need to be positive and "
            "monotonically increasing.");
    }
}

void PlasmaPhase::checkElectronEnergyDistribution() const
{
    Eigen::ArrayXd h = m_electronEnergyLevels.tail(m_nPoints - 1) -
                       m_electronEnergyLevels.head(m_nPoints - 1);
    if ((m_electronEnergyDist < 0.0).any()) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
            "Values of electron energy distribution cannot be negative.");
    }
    if (m_electronEnergyDist[m_nPoints - 1] > 0.01) {
        warn_user("PlasmaPhase::checkElectronEnergyDistribution",
        "The value of the last element of electron energy distribution exceed 0.01. "
        "This indicates that the value of electron energy level is not high enough "
        "to contain the isotropic distribution at mean electron energy of "
        "{} eV", meanElectronEnergy());
    }
}

void PlasmaPhase::setDiscretizedElectronEnergyDist(const double* levels,
                                                const double* dist,
                                                size_t length)
{
    m_distributionType = "discretized";
    m_nPoints = length;
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::ArrayXd>(levels, length);
    m_electronEnergyDist =
        Eigen::Map<const Eigen::ArrayXd>(dist, length);
    checkElectronEnergyLevels();
    if (m_do_normalizeElectronEnergyDist) {
        normalizeElectronEnergyDistribution();
    }
    checkElectronEnergyDistribution();
    updateElectronEnergyDistDifference();
    updateElectronTemperatureFromEnergyDist();
    electronEnergyLevelChanged();
    electronEnergyDistributionChanged();
}

void PlasmaPhase::updateElectronTemperatureFromEnergyDist()
{
    // calculate mean electron energy and electron temperature
    Eigen::ArrayXd eps52 = m_electronEnergyLevels.pow(5./2.);
    double epsilon_m = 2.0 / 5.0 * numericalQuadrature(m_quadratureMethod,
                                                       m_electronEnergyDist, eps52);
    m_electronTemp = 2.0 / 3.0 * epsilon_m * ElectronCharge / Boltzmann;
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
    vector<double> levels(m_nPoints);
    Eigen::Map<Eigen::ArrayXd>(levels.data(), m_nPoints) = m_electronEnergyLevels;
    eedf["energy-levels"] = levels;
    if (m_distributionType == "isotropic") {
        eedf["shape-factor"] = m_isotropicShapeFactor;
        eedf["mean-electron-energy"].setQuantity(meanElectronEnergy(), "eV");
    } else if (m_distributionType == "discretized") {
        vector<double> dist(m_nPoints);
        Eigen::Map<Eigen::ArrayXd>(dist.data(), m_nPoints) = m_electronEnergyDist;
        eedf["distribution"] = dist;
        eedf["normalize"] = m_do_normalizeElectronEnergyDist;
    }
    phaseNode["electron-energy-distribution"] = std::move(eedf);
}

void PlasmaPhase::setParameters(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    IdealGasPhase::setParameters(phaseNode, rootNode);
    m_root = rootNode;
    if (phaseNode.hasKey("electron-energy-distribution")) {
        const AnyMap eedf = phaseNode["electron-energy-distribution"].as<AnyMap>();
        m_distributionType = eedf["type"].asString();
        if (m_distributionType == "isotropic") {
            if (eedf.hasKey("shape-factor")) {
                setIsotropicShapeFactor(eedf["shape-factor"].asDouble());
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "isotropic type requires shape-factor key.");
            }
            if (eedf.hasKey("mean-electron-energy")) {
                double energy = eedf.convert("mean-electron-energy", "eV");
                setMeanElectronEnergy(energy);
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "isotropic type requires electron-temperature key.");
            }
            if (eedf.hasKey("energy-levels")) {
                setElectronEnergyLevels(eedf["energy-levels"].asVector<double>().data(),
                                        eedf["energy-levels"].asVector<double>().size());
            }
            setIsotropicElectronEnergyDistribution();
        } else if (m_distributionType == "discretized") {
            if (!eedf.hasKey("energy-levels")) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cannot find key energy-levels.");
            }
            if (!eedf.hasKey("distribution")) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cannot find key distribution.");
            }
            if (eedf.hasKey("normalize")) {
                enableNormalizeElectronEnergyDist(eedf["normalize"].asBool());
            }
            setDiscretizedElectronEnergyDist(eedf["energy-levels"].asVector<double>().data(),
                                             eedf["distribution"].asVector<double>().data(),
                                             eedf["energy-levels"].asVector<double>().size());
        }
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

void PlasmaPhase::setSolution(std::weak_ptr<Solution> soln) {
    ThermoPhase::setSolution(soln);
    // register callback function to be executed
    // when the thermo or kinetics object changed
    if (shared_ptr<Solution> soln = m_soln.lock()) {
        soln->registerChangedCallback(this, [&]() {
            setCollisions();
        });
    }
}

void PlasmaPhase::setCollisions()
{
    m_collisions.clear();
    m_collisionRates.clear();
    m_targetSpeciesIndices.clear();

    if (shared_ptr<Solution> soln = m_soln.lock()) {
        shared_ptr<Kinetics> kin = soln->kinetics();
        if (!kin || kin->nReactions() == 0) {
            return;
        }

        // Collect electron-collision-plasma reactions
        for (size_t i = 0; i < kin->nReactions(); i++) {
            std::shared_ptr<Reaction> R = kin->reaction(i);
            if (R->rate()->type() != "electron-collision-plasma") {
                continue;
            }
            m_collisions.push_back(R);
        }

        // Resize based on the size of m_collisions
        m_collisionRates.resize(nCollisions());
        m_targetSpeciesIndices.resize(nCollisions());

        for (size_t i = 0; i < nCollisions(); i++) {
            std::shared_ptr<Reaction> R = m_collisions[i];
            m_collisionRates[i] = boost::polymorphic_pointer_downcast
                                    <ElectronCollisionPlasmaRate>(R->rate());

            // Identify target species for electron-collision reactions
            for (const auto& [name, stoich] : R->reactants) {
                // Reactants are expected to be electrons and the target species
                if (name != electronSpeciesName()) {
                    m_targetSpeciesIndices[i] = speciesIndex(name);
                    break;
                }
            }
        }
        updateInterpolatedCrossSections();
    }
}

void PlasmaPhase::updateInterpolatedCrossSections()
{
    for (shared_ptr<ElectronCollisionPlasmaRate> rate : m_collisionRates) {
        for (size_t i = 0; i < nElectronEnergyLevels(); i++) {
            double level = m_electronEnergyLevels[i];
            m_interp_cs[i] = linearInterp(level,
                rate->energyLevels(), rate->crossSections());
        }
        // Set the interpolated cross section
        rate->setCrossSectionInterpolated(m_interp_cs);
    }
    interpolatedCrossSectionsChanged();
}

void PlasmaPhase::updateElectronEnergyDistDifference()
{
    m_electronEnergyDistDiff.resize(nElectronEnergyLevels());
    // Forward difference for the first point
    m_electronEnergyDistDiff[0] =
        (m_electronEnergyDist[1] - m_electronEnergyDist[0]) /
        (m_electronEnergyLevels[1] - m_electronEnergyLevels[0]);

    // Central difference for the middle points
    for (size_t i = 1; i < m_nPoints - 1; i++) {
        double h1 = m_electronEnergyLevels[i+1] - m_electronEnergyLevels[i];
        double h0 = m_electronEnergyLevels[i] - m_electronEnergyLevels[i-1];
        m_electronEnergyDistDiff[i] = (h0 * h0 * m_electronEnergyDist[i+1] +
                (h1 * h1 - h0 * h0) * m_electronEnergyDist[i] -
                h1 * h1 * m_electronEnergyDist[i-1]) /
                (h1 * h0) / (h1 + h0);
    }

    // Backward difference for the last point
    m_electronEnergyDistDiff[m_nPoints-1] =
        (m_electronEnergyDist[m_nPoints-1] -
        m_electronEnergyDist[m_nPoints-2]) /
        (m_electronEnergyLevels[m_nPoints-1] -
        m_electronEnergyLevels[m_nPoints-2]);
}

void PlasmaPhase::updateElasticElectronEnergyLossCoefficients()
{
    static const int cacheId = m_cache.getId();
    CachedScalar last = m_cache.getScalar(cacheId);

    // both change of distribution and cross section cause
    // the elastic electron energy loss coefficients change
    int statenum = m_csNum + m_distNum;

    if (last.validate(temperature(), statenum)) return;

    m_elasticElectronEnergyLossCoefficients.resize(nCollisions());

    // @todo exclude attachment collisions
    for (size_t i = 0; i < nCollisions(); i++) {
        size_t k = m_targetSpeciesIndices[i];

        // Map cross sections to Eigen::ArrayXd
        auto cs_array = Eigen::Map<const Eigen::ArrayXd>(
            m_collisionRates[i]->crossSectionInterpolated().data(),
            m_collisionRates[i]->crossSectionInterpolated().size()
        );

        // Mass ratio calculation
        double mass_ratio = ElectronMass / molecularWeight(k) * Avogadro;

        // Calculate the rate using Simpson's rule
        Eigen::ArrayXd f0_plus = m_electronEnergyDist + Boltzmann * temperature() /
                                 ElectronCharge * m_electronEnergyDistDiff;
        m_elasticElectronEnergyLossCoefficients[i] = 2.0 * mass_ratio * gamma *
            simpson(1.0 / 3.0 * f0_plus.cwiseProduct(cs_array)
                    , m_electronEnergyLevels.pow(3.0));
    }
}

double PlasmaPhase::elasticPowerLoss()
{
    updateElasticElectronEnergyLossCoefficients();

    // The elastic power loss includes the contributions from inelastic
    // collisions (inelastic recoil effects).
    double rate = 0.0;
    for (size_t i = 0; i < nCollisions(); i++) {
        rate += concentration(m_targetSpeciesIndices[i]) *
                m_elasticElectronEnergyLossCoefficients[i];
    }

    return Avogadro * Avogadro * ElectronCharge *
        concentration(m_electronSpeciesIndex) * rate;
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

double PlasmaPhase::enthalpy_mole() const {
    double value = IdealGasPhase::enthalpy_mole();
    value += GasConstant * (electronTemperature() - temperature()) *
             moleFraction(m_electronSpeciesIndex) *
             m_h0_RT[m_electronSpeciesIndex];
    return value;
}

void PlasmaPhase::getGibbs_ref(double* g) const
{
    IdealGasPhase::getGibbs_ref(g);
    g[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getStandardVolumes_ref(double* vol) const
{
    IdealGasPhase::getStandardVolumes_ref(vol);
    vol[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEnthalpies(double* hbar) const
{
    IdealGasPhase::getPartialMolarEnthalpies(hbar);
    hbar[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEntropies(double* sbar) const
{
    IdealGasPhase::getPartialMolarEntropies(sbar);
    double logp = log(pressure());
    double logpe = log(electronPressure());
    sbar[m_electronSpeciesIndex] += GasConstant * (logp - logpe);
}

void PlasmaPhase::getPartialMolarIntEnergies(double* ubar) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
    size_t k = m_electronSpeciesIndex;
    ubar[k] = RTe() * (_h[k] - 1.0);
}

void PlasmaPhase::getChemPotentials(double* mu) const
{
    IdealGasPhase::getChemPotentials(mu);
    size_t k = m_electronSpeciesIndex;
    double xx = std::max(SmallNumber, moleFraction(k));
    mu[k] += (RTe() - RT()) * log(xx);
}

void PlasmaPhase::getStandardChemPotentials(double* muStar) const
{
    IdealGasPhase::getStandardChemPotentials(muStar);
    size_t k = m_electronSpeciesIndex;
    muStar[k] -= log(pressure() / refPressure()) * RT();
    muStar[k] += log(electronPressure() / refPressure()) * RTe();
}

void PlasmaPhase::getEntropy_R(double* sr) const
{
    const vector<double>& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            sr[k] -= tmp;
        } else {
            sr[k] -= log(electronPressure() / refPressure());
        }
    }
}

void PlasmaPhase::getGibbs_RT(double* grt) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            grt[k] += tmp;
        } else {
            grt[k] += log(electronPressure() / refPressure());
        }
    }
}

}
