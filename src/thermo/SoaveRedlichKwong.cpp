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

SoaveRedlichKwong::SoaveRedlichKwong(const std::string& infile, const std::string& id_)  :
    m_b(0.0),
    m_a(0.0),
    m_aAlpha_mix(0.0),
    m_NSolns(0),
    m_dpdV(0.0),
    m_dpdT(0.0)
{
    std::fill_n(m_Vroot, 3, 0.0);
    initThermoFile(infile, id_);
}

void SoaveRedlichKwong::setSpeciesCoeffs(const std::string& species, double a, double b,
                                         double w)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
        throw CanteraError("SoaveRedlichKwong::setSpeciesCoeffs",
            "Unknown species '{}'.", species);
    }

    // Calculate value of kappa (independent of temperature)
    // w is an acentric factor of species
    m_kappa[k] = 0.480 + 1.574*w - 0.176*w*w;
    m_acentric[k] = w; // store the original acentric factor to enable serialization

    // Calculate alpha (temperature dependent interaction parameter)
    double critTemp = speciesCritTemperature(a, b);
    double sqt_T_r = sqrt(temperature() / critTemp);
    double sqt_alpha = 1 + m_kappa[k] * (1 - sqt_T_r);
    m_alpha[k] = sqt_alpha*sqt_alpha;

    m_a_coeffs(k,k) = a;
    double aAlpha_k = a*m_alpha[k];
    m_aAlpha_binary(k,k) = aAlpha_k;

    // standard mixing rule for cross-species interaction term
    for (size_t j = 0; j < m_kk; j++) {
        if (k == j) {
            continue;
        }
        double a0kj = sqrt(m_a_coeffs(j,j) * a);
        double aAlpha_j = a*m_alpha[j];
        double a_Alpha = sqrt(aAlpha_j*aAlpha_k);
        if (m_a_coeffs(j, k) == 0) {
            m_a_coeffs(j, k) = a0kj;
            m_aAlpha_binary(j, k) = a_Alpha;
            m_a_coeffs(k, j) = a0kj;
            m_aAlpha_binary(k, j) = a_Alpha;
        }
    }
    m_b_coeffs[k] = b;
}

void SoaveRedlichKwong::setBinaryCoeffs(const std::string& species_i,
        const std::string& species_j, double a0)
{
    size_t ki = speciesIndex(species_i);
    if (ki == npos) {
        throw CanteraError("SoaveRedlichKwong::setBinaryCoeffs",
            "Unknown species '{}'.", species_i);
    }
    size_t kj = speciesIndex(species_j);
    if (kj == npos) {
        throw CanteraError("SoaveRedlichKwong::setBinaryCoeffs",
            "Unknown species '{}'.", species_j);
    }

    m_a_coeffs(ki, kj) = m_a_coeffs(kj, ki) = a0;
    m_binaryParameters[species_i][species_j] = a0;
    m_binaryParameters[species_j][species_i] = a0;
    // Calculate alpha_ij
    double alpha_ij = m_alpha[ki] * m_alpha[kj];
    m_aAlpha_binary(ki, kj) = m_aAlpha_binary(kj, ki) = a0*alpha_ij;
}

// ------------Molar Thermodynamic Properties -------------------------

double SoaveRedlichKwong::cp_mole() const
{
    return 1.0;
}

double SoaveRedlichKwong::cv_mole() const
{
    _updateReferenceStateThermo();
    double T = temperature();
    calculatePressureDerivatives();
    return (cp_mole() + T * m_dpdT * m_dpdT / m_dpdV);
}

double SoaveRedlichKwong::pressure() const
{
    _updateReferenceStateThermo();
    // Get a copy of the private variables stored in the State object
    double T = temperature();
    double mv = molarVolume();
    return GasConstant * T / (mv - m_b) - m_aAlpha_mix / (mv * (mv + m_b));
}

double SoaveRedlichKwong::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void SoaveRedlichKwong::getActivityCoefficients(double* ac) const
{
    
}

// ---- Partial Molar Properties of the Solution -----------------

void SoaveRedlichKwong::getChemPotentials(double* mu) const
{

}

void SoaveRedlichKwong::getPartialMolarEnthalpies(double* hbar) const
{

}

void SoaveRedlichKwong::getPartialMolarEntropies(double* sbar) const
{

}

void SoaveRedlichKwong::getPartialMolarIntEnergies(double* ubar) const
{

}

void SoaveRedlichKwong::getPartialMolarCp(double* cpbar) const
{

}

void SoaveRedlichKwong::getPartialMolarVolumes(double* vbar) const
{

}

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

double SoaveRedlichKwong::sresid() const
{ // Replace
    double molarV = molarVolume();
    double hh = m_b / molarV;
    double zz = z();
    double alpha_1 = daAlpha_dT();
    double vpb = molarV + (1.0 + Sqrt2) * m_b;
    double vmb = molarV + (1.0 - Sqrt2) * m_b;
    double fac = alpha_1 / (2.0 * Sqrt2 * m_b);
    double sresid_mol_R = log(zz*(1.0 - hh)) + fac * log(vpb / vmb) / GasConstant;
    return GasConstant * sresid_mol_R;
}

double SoaveRedlichKwong::hresid() const
{ // Replace
    double molarV = molarVolume();
    double zz = z();
    double aAlpha_1 = daAlpha_dT();
    double T = temperature();
    double vpb = molarV + (1 + Sqrt2) * m_b;
    double vmb = molarV + (1 - Sqrt2) * m_b;
    double fac = 1 / (2.0 * Sqrt2 * m_b);
    return GasConstant * T * (zz - 1.0)
        + fac * log(vpb / vmb) * (T * aAlpha_1 - m_aAlpha_mix);
}

double SoaveRedlichKwong::liquidVolEst(double T, double& presGuess) const
{
    double v = m_b * 1.1;
    double atmp, btmp, aAlphatmp;
    calculateAB(atmp, btmp, aAlphatmp);
    double pres = std::max(psatEst(T), presGuess);
    double Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
        int nsol = solveCubic(T, pres, atmp, btmp, aAlphatmp, Vroot);
        if (nsol == 1 || nsol == 2) {
            double pc = critPressure();
            if (pres > pc) {
                foundLiq = true;
            }
            pres *= 1.04;
        } else {
            foundLiq = true;
        }
    }

    if (foundLiq) {
        v = Vroot[0];
        presGuess = pres;
    } else {
        v = -1.0;
    }
    return v;
}

double SoaveRedlichKwong::densityCalc(double T, double presPa, int phaseRequested,
                                 double rhoGuess)
{
    // It's necessary to set the temperature so that m_aAlpha_mix is set correctly.
    setTemperature(T);
    double tcrit = critTemperature();
    double mmw = meanMolecularWeight();
    if (rhoGuess == -1.0) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            double lqvol = liquidVolEst(T, presPa);
            rhoGuess = mmw / lqvol;
        }
    } else {
        // Assume the Gas phase initial guess, if nothing is specified to the routine
        rhoGuess = presPa * mmw / (GasConstant * T);
    }

    double volGuess = mmw / rhoGuess;
    m_NSolns = solveCubic(T, presPa, m_a, m_b, m_aAlpha_mix, m_Vroot);

    double molarVolLast = m_Vroot[0];
    if (m_NSolns >= 2) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            molarVolLast = m_Vroot[0];
        } else if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = m_Vroot[2];
        } else {
            if (volGuess > m_Vroot[1]) {
                molarVolLast = m_Vroot[2];
            } else {
                molarVolLast = m_Vroot[0];
            }
        }
    } else if (m_NSolns == 1) {
        if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT
            || phaseRequested == FLUID_UNDEFINED)
        {
            molarVolLast = m_Vroot[0];
        } else {
            return -2.0;
        }
    } else if (m_NSolns == -1) {
        if (phaseRequested >= FLUID_LIQUID_0 || phaseRequested == FLUID_UNDEFINED
            || phaseRequested == FLUID_SUPERCRIT)
        {
            molarVolLast = m_Vroot[0];
        } else if (T > tcrit) {
            molarVolLast = m_Vroot[0];
        } else {
            return -2.0;
        }
    } else {
        molarVolLast = m_Vroot[0];
        return -1.0;
    }
    return mmw / molarVolLast;
}

double SoaveRedlichKwong::densSpinodalLiquid() const
{
    double Vroot[3];
    double T = temperature();
    int nsol = solveCubic(T, pressure(), m_a, m_b, m_aAlpha_mix, Vroot);
    if (nsol != 3) {
        return critDensity();
    }

    auto resid = [this, T](double v) {
        double pp;
        return dpdVCalc(T, v, pp);
    };

    boost::uintmax_t maxiter = 100;
    std::pair<double, double> vv = bmt::toms748_solve(
        resid, Vroot[0], Vroot[1], bmt::eps_tolerance<double>(48), maxiter);

    double mmw = meanMolecularWeight();
    return mmw / (0.5 * (vv.first + vv.second));
}

double SoaveRedlichKwong::densSpinodalGas() const
{
    double Vroot[3];
    double T = temperature();
    int nsol = solveCubic(T, pressure(), m_a, m_b, m_aAlpha_mix, Vroot);
    if (nsol != 3) {
        return critDensity();
    }

    auto resid = [this, T](double v) {
        double pp;
        return dpdVCalc(T, v, pp);
    };

    boost::uintmax_t maxiter = 100;
    std::pair<double, double> vv = bmt::toms748_solve(
        resid, Vroot[1], Vroot[2], bmt::eps_tolerance<double>(48), maxiter);

    double mmw = meanMolecularWeight();
    return mmw / (0.5 * (vv.first + vv.second));
}

double SoaveRedlichKwong::dpdVCalc(double T, double molarVol, double& presCalc) const
{
    double denom = molarVol * (molarVol + m_b);
    double vmb = molarVol - m_b;
    return -GasConstant * T / (vmb * vmb) + m_aAlpha_mix * (2 * molarVol + m_b) / (denom * denom);
}

double SoaveRedlichKwong::isothermalCompressibility() const
{
    calculatePressureDerivatives();
    return -1 / (molarVolume() * m_dpdV);
}

double SoaveRedlichKwong::thermalExpansionCoeff() const
{
    calculatePressureDerivatives();
    return -m_dpdT / (molarVolume() * m_dpdV);
}

void SoaveRedlichKwong::calculatePressureDerivatives() const
{
    double T = temperature();
    double mv = molarVolume();
    double pres;

    m_dpdV = dpdVCalc(T, mv, pres);
    m_dpdT = GasConstant / (mv - m_b) - daAlpha_dT() / (mv * (mv + m_b));
}

void SoaveRedlichKwong::updateMixingExpressions()
{
    // Same as PengRobinson::updateMixingExpressions()
    double temp = temperature();

    // Update individual alpha
    for (size_t j = 0; j < m_kk; j++) {
        double critTemp_j = speciesCritTemperature(m_a_coeffs(j,j), m_b_coeffs[j]);
        double sqt_alpha = 1 + m_kappa[j] * (1 - sqrt(temp / critTemp_j));
        m_alpha[j] = sqt_alpha*sqt_alpha;
    }

    // Update aAlpha_i, j
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            m_aAlpha_binary(i, j) = sqrt(m_alpha[i] * m_alpha[j]) * m_a_coeffs(i,j);
        }
    }
    calculateAB(m_a,m_b,m_aAlpha_mix);
}

void SoaveRedlichKwong::calculateAB(double& aCalc, double& bCalc, double& aAlphaCalc) const
{
    // Same as PengRobinson::calculateAB()
    bCalc = 0.0;
    aCalc = 0.0;
    aAlphaCalc = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        bCalc += moleFractions_[i] * m_b_coeffs[i];
        for (size_t j = 0; j < m_kk; j++) {
            aCalc += m_a_coeffs(i, j) * moleFractions_[i] * moleFractions_[j];
            aAlphaCalc += m_aAlpha_binary(i, j) * moleFractions_[i] * moleFractions_[j];
        }
    }
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

void SoaveRedlichKwong::calcCriticalConditions(double& pc, double& tc, double& vc) const
{
    if (m_b <= 0.0) {
        tc = 1000000.;
        pc = 1.0E13;
        vc = omega_vc * GasConstant * tc / pc;
        return;
    }
    if (m_a <= 0.0) {
        tc = 0.0;
        pc = 0.0;
        vc = 2.0 * m_b;
        return;
    }
    tc = m_a * omega_b / (m_b * omega_a * GasConstant);
    pc = omega_b * GasConstant * tc / m_b;
    vc = omega_vc * GasConstant * tc / pc;
}

int SoaveRedlichKwong::solveCubic(double T, double pres, double a, double b, double aAlpha,
                                  double Vroot[3]) const
{
    double an = 1.0;
    double bn = - GasConstant * T / pres;
    double cn = (aAlpha - b * GasConstant * T) / pres - b * b;
    double dn = -aAlpha * b / pres;

    double tc = a * omega_b / (b * omega_a * GasConstant);
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    return MixtureFugacityTP::solveCubic(T, pres, a, b, aAlpha, Vroot,
                                         an, bn, cn, dn, tc, vc);
}

}
