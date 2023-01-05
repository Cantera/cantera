//! @file SoaveRedlichKwong.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SoaveRedlichKwong.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"

#include <boost/math/tools/roots.hpp>

namespace bmt = boost::math::tools;

namespace Cantera
{
const double SoaveRedlichKwong::omega_a = 4.27480233540E-01;
const double SoaveRedlichKwong::omega_b = 8.66403499650E-02;
const double SoaveRedlichKwong::omega_vc = 3.33333333333333E-01;


double SoaveRedlichKwong::speciesCritTemperature(double a, double b) const
{
    if (b <= 0.0) {
        return 1000000.;
    } else if (a <= 0.0) {
        return 0.0;
    } else {
        return a * omega_b / (b * omega_a * GasConstant);
    }
}


bool SoaveRedlichKwong::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        m_a_coeffs.resize(m_kk, m_kk, 0.0);
        m_b_coeffs.push_back(0.0);
        m_aAlpha_binary.resize(m_kk, m_kk, 0.0);
        m_kappa.push_back(0.0);
        m_acentric.push_back(0.0);
        m_alpha.push_back(0.0);
        m_dalphadT.push_back(0.0);
        m_d2alphadT2.push_back(0.0);
        m_pp.push_back(0.0);
        m_partialMolarVolumes.push_back(0.0);
        m_dpdni.push_back(0.0);
        m_coeffSource.push_back(CoeffSource::EoS);
    }
    return added;
}


void SoaveRedlichKwong::initThermo()
{
    // Contents of 'critical-properties.yaml', loaded later if needed
    AnyMap critPropsDb;
    std::unordered_map<std::string, AnyMap*> dbSpecies;

    for (auto& item : m_species) {
        auto& data = item.second->input;
        size_t k = speciesIndex(item.first);
        if (m_a_coeffs(k, k) != 0.0) {
            continue;
        }
        bool foundCoeffs = false;
        if (data.hasKey("equation-of-state") &&
            data["equation-of-state"].hasMapWhere("model", "Soave-Redlich-Kwong"))
        {
            // Read a and b coefficients and acentric factor w_ac from species input
            // information, specified in a YAML input file.
            auto eos = data["equation-of-state"].getMapWhere(
                "model", "Soave-Redlich-Kwong");
            if (eos.hasKey("a") && eos.hasKey("b") && eos.hasKey("acentric-factor")) {
                double a0 = eos.convert("a", "Pa*m^6/kmol^2");
                double b = eos.convert("b", "m^3/kmol");
                // unitless acentric factor:
                double w = eos["acentric-factor"].asDouble();
                setSpeciesCoeffs(item.first, a0, b, w);
                foundCoeffs = true;
            }

            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& item2 : binary_a) {
                    double a0 = units.convert(item2.second, "Pa*m^6/kmol^2");
                    setBinaryCoeffs(item.first, item2.first, a0);
                }
            }
            if (foundCoeffs) {
                m_coeffSource[k] = CoeffSource::EoS;
                continue;
            }
        }

        // Coefficients have not been populated from model-specific input
        double Tc = NAN, Pc = NAN, omega_ac = NAN;
        if (data.hasKey("critical-parameters")) {
            // Use critical state information stored in the species entry to
            // calculate a, b, and the acentric factor.
            auto& critProps = data["critical-parameters"].as<AnyMap>();
            Tc = critProps.convert("critical-temperature", "K");
            Pc = critProps.convert("critical-pressure", "Pa");
            omega_ac = critProps["acentric-factor"].asDouble();
            m_coeffSource[k] = CoeffSource::CritProps;
        } else {
            // Search 'crit-properties.yaml' to find Tc and Pc. Load data if needed.
            if (critPropsDb.empty()) {
                critPropsDb = AnyMap::fromYamlFile("critical-properties.yaml");
                dbSpecies = critPropsDb["species"].asMap("name");
            }

            // All names in critical-properties.yaml are upper case
            auto ucName = boost::algorithm::to_upper_copy(item.first);
            if (dbSpecies.count(ucName)) {
                auto& spec = *dbSpecies.at(ucName);
                auto& critProps = spec["critical-parameters"].as<AnyMap>();
                Tc = critProps.convert("critical-temperature", "K");
                Pc = critProps.convert("critical-pressure", "Pa");
                omega_ac = critProps["acentric-factor"].asDouble();
                m_coeffSource[k] = CoeffSource::Database;
            }
        }

        // Check if critical properties were found in either location
        if (!isnan(Tc)) {
            double a = omega_a * std::pow(GasConstant * Tc, 2) / Pc;
            double b = omega_b * GasConstant * Tc / Pc;
            setSpeciesCoeffs(item.first, a, b, omega_ac);
        } else {
            throw InputFileError("SoaveRedlichKwong::initThermo", data,
            "No Soave-Redlich-Kwong model parameters or critical properties found for "
            "species '{}'", item.first);
        }
    }
}


void SoaveRedlichKwong::getSpeciesParameters(const std::string& name,
                                        AnyMap& speciesNode) const
{
    MixtureFugacityTP::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name);
    checkSpeciesIndex(k);

    // Pure species parameters
    if (m_coeffSource[k] == CoeffSource::EoS) {
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Soave-Redlich-Kwong", true);
        eosNode["a"].setQuantity(m_a_coeffs(k, k), "Pa*m^6/kmol^2");
        eosNode["b"].setQuantity(m_b_coeffs[k], "m^3/kmol");
        eosNode["acentric-factor"] = m_acentric[k];
    } else if (m_coeffSource[k] == CoeffSource::CritProps) {
        auto& critProps = speciesNode["critical-parameters"];
        double Tc = speciesCritTemperature(m_a_coeffs(k, k), m_b_coeffs[k]);
        double Pc = omega_b * GasConstant * Tc / m_b_coeffs[k];
        critProps["critical-temperature"].setQuantity(Tc, "K");
        critProps["critical-pressure"].setQuantity(Pc, "Pa");
        critProps["acentric-factor"] = m_acentric[k];
    }
    // Nothing to do in the case where the parameters are from the database

    if (m_binaryParameters.count(name)) {
        // Include binary parameters regardless of where the pure species parameters
        // were found
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Soave-Redlich-Kwong", true);
        AnyMap bin_a;
        for (const auto& item : m_binaryParameters.at(name)) {
            bin_a[item.first].setQuantity(item.second, "Pa*m^6/kmol^2");
        }
        eosNode["binary-a"] = std::move(bin_a);
    }
}

double SoaveRedlichKwong::dpdVCalc(double T, double molarVol, double& presCalc) const
{
    double denom = molarVol * (molarVol + m_b);
    double vmb = molarVol - m_b;
    return -GasConstant * T / (vmb * vmb) + m_aAlpha_mix * (2 * molarVol + m_b) / (denom * denom);
}

void SoaveRedlichKwong::calculatePressureDerivatives() const
{
    double T = temperature();
    double mv = molarVolume();
    double pres;

    m_dpdV = dpdVCalc(T, mv, pres);
    m_dpdT = GasConstant / (mv - m_b) - daAlpha_dT() / (mv * (mv + m_b));
}

double SoaveRedlichKwong::daAlpha_dT() const
{
    // Same as PengRobinson::daAlpha_dT()
    double daAlphadT = 0.0, k, Tc, sqtTr, coeff1, coeff2;
    for (size_t i = 0; i < m_kk; i++) {
        // Calculate first derivative of alpha for individual species
        Tc = speciesCritTemperature(m_a_coeffs(i,i), m_b_coeffs[i]);
        sqtTr = sqrt(temperature() / Tc);
        coeff1 = 1 / (Tc*sqtTr);
        coeff2 = sqtTr - 1;
        k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
    }
    // Calculate mixture derivative
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            daAlphadT += moleFractions_[i] * moleFractions_[j] * 0.5
                         * m_aAlpha_binary(i, j)
                         * (m_dalphadT[i] / m_alpha[i] + m_dalphadT[j] / m_alpha[j]);
        }
    }
    return daAlphadT;
}

double SoaveRedlichKwong::d2aAlpha_dT2() const
{
    // Same as PengRobinson::d2aAlpha_dT2()
    for (size_t i = 0; i < m_kk; i++) {
        double Tcrit_i = speciesCritTemperature(m_a_coeffs(i, i), m_b_coeffs[i]);
        double sqt_Tr = sqrt(temperature() / Tcrit_i);
        double coeff1 = 1 / (Tcrit_i*sqt_Tr);
        double coeff2 = sqt_Tr - 1;
        // Calculate first and second derivatives of alpha for individual species
        double k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
        m_d2alphadT2[i] = (k*k + k) * coeff1 / (2*sqt_Tr*sqt_Tr*Tcrit_i);
    }

    // Calculate mixture derivative
    double d2aAlphadT2 = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        double alphai = m_alpha[i];
        for (size_t j = 0; j < m_kk; j++) {
            double alphaj = m_alpha[j];
            double alphaij = alphai * alphaj;
            double term1 = m_d2alphadT2[i] / alphai + m_d2alphadT2[j] / alphaj;
            double term2 = 2 * m_dalphadT[i] * m_dalphadT[j] / alphaij;
            double term3 = m_dalphadT[i] / alphai + m_dalphadT[j] / alphaj;
            d2aAlphadT2 += 0.5 * moleFractions_[i] * moleFractions_[j]
                           * m_aAlpha_binary(i, j)
                           * (term1 + term2 - 0.5 * term3 * term3);
        }
    }
    return d2aAlphadT2;
}


int SoaveRedlichKwong::solveCubic(double T, double pres, double a, double b, double aAlpha,
                                  double Vroot[3]) const
{
    double an = 1.0;
    double bn = - GasConstant * T / pres;
    double cn = (aAlpha - b * GasConstant * T) / pres - b * b;
    double dn = aAlpha * b / pres;

    double tc = a * omega_b / (b * omega_a * GasConstant);
    double pc = omega_b * R * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    return MixtureFugacityTP::solveCubic(T, pres, a, b, aAlpha, Vroot,
                                         an, bn, cn, dn, tc, vc);
}

}
