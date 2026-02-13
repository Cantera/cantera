//! @file RedlichKwongMFTP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

#include <boost/algorithm/string.hpp>
#include <algorithm>

namespace Cantera
{

const double RedlichKwongMFTP::omega_a = 4.27480233540E-01;
const double RedlichKwongMFTP::omega_b = 8.66403499650E-02;
const double RedlichKwongMFTP::omega_vc = 3.33333333333333E-01;

RedlichKwongMFTP::RedlichKwongMFTP(const string& infile, const string& id_)
{
    initThermoFile(infile, id_);
}

void RedlichKwongMFTP::setSpeciesCoeffs(const string& species,
                                        double a0, double a1, double b)
{
    size_t k = speciesIndex(species, true);

    if (a1 != 0.0) {
        m_formTempParam = 1; // expression is temperature-dependent
    }

    size_t counter = k + m_kk * k;
    a_coeff_vec(0, counter) = a0;
    a_coeff_vec(1, counter) = a1;

    // standard mixing rule for cross-species interaction term
    for (size_t j = 0; j < m_kk; j++) {
        if (k == j) {
            continue;
        }

        // a_coeff_vec(0) is initialized to NaN to mark uninitialized species
        if (isnan(a_coeff_vec(0, j + m_kk * j))) {
            // The diagonal element of the jth species has not yet been defined.
            continue;
        } else if (isnan(a_coeff_vec(0, j + m_kk * k))) {
            // Only use the mixing rules if the off-diagonal element has not already been defined by a
            // user-specified crossFluidParameters entry:
            double a0kj = sqrt(a_coeff_vec(0, j + m_kk * j) * a0);
            double a1kj = sqrt(a_coeff_vec(1, j + m_kk * j) * a1);
            a_coeff_vec(0, j + m_kk * k) = a0kj;
            a_coeff_vec(1, j + m_kk * k) = a1kj;
            a_coeff_vec(0, k + m_kk * j) = a0kj;
            a_coeff_vec(1, k + m_kk * j) = a1kj;
        }
    }
    a_coeff_vec.getRow(0, a_vec_Curr_);
    b_vec_Curr_[k] = b;
}

void RedlichKwongMFTP::setBinaryCoeffs(const string& species_i, const string& species_j,
                                       double a0, double a1)
{
    size_t ki = speciesIndex(species_i, true);
    size_t kj = speciesIndex(species_j, true);

    if (a1 != 0.0) {
        m_formTempParam = 1; // expression is temperature-dependent
    }

    m_binaryParameters[species_i][species_j] = {a0, a1};
    m_binaryParameters[species_j][species_i] = {a0, a1};
    size_t counter1 = ki + m_kk * kj;
    size_t counter2 = kj + m_kk * ki;
    a_coeff_vec(0, counter1) = a_coeff_vec(0, counter2) = a0;
    a_coeff_vec(1, counter1) = a_coeff_vec(1, counter2) = a1;
    a_vec_Curr_[counter1] = a_vec_Curr_[counter2] = a0;
}

// ------------Molar Thermodynamic Properties -------------------------

double RedlichKwongMFTP::cp_mole() const
{
    _updateReferenceStateThermo();
    double TKelvin = temperature();
    double sqt = sqrt(TKelvin);
    double mv = molarVolume();
    double vpb = mv + m_b_current;
    pressureDerivatives();
    double cpref = GasConstant * mean_X(m_cp0_R);
    double dadt = da_dt();
    double fac = TKelvin * dadt - 3.0 * m_a_current / 2.0;
    double dHdT_V = (cpref + mv * dpdT_ - GasConstant - 1.0 / (2.0 * m_b_current * TKelvin * sqt) * log(vpb/mv) * fac
                         +1.0/(m_b_current * sqt) * log(vpb/mv) * (-0.5 * dadt));
    return dHdT_V - (mv + TKelvin * dpdT_ / dpdV_) * dpdT_;
}

double RedlichKwongMFTP::cv_mole() const
{
    _updateReferenceStateThermo();
    double TKelvin = temperature();
    double sqt = sqrt(TKelvin);
    double mv = molarVolume();
    double vpb = mv + m_b_current;
    double cvref = GasConstant * (mean_X(m_cp0_R) - 1.0);
    double dadt = da_dt();
    double fac = TKelvin * dadt - 3.0 * m_a_current / 2.0;
    return (cvref - 1.0/(2.0 * m_b_current * TKelvin * sqt) * log(vpb/mv)*fac
            +1.0/(m_b_current * sqt) * log(vpb/mv)*(-0.5*dadt));
}

double RedlichKwongMFTP::pressure() const
{
    _updateReferenceStateThermo();

    //  Get a copy of the private variables stored in the State object
    double T = temperature();
    double molarV = meanMolecularWeight() / density();
    double pp = GasConstant * T/(molarV - m_b_current) - m_a_current/(sqrt(T) * molarV * (molarV + m_b_current));
    return pp;
}

double RedlichKwongMFTP::standardConcentration(size_t k) const
{
    getStandardVolumes(m_workS);
    return 1.0 / m_workS[k];
}

void RedlichKwongMFTP::getActivityCoefficients(span<double> ac) const
{
    checkArraySize("RedlichKwongMFTP::getActivityCoefficients", ac.size(), m_kk);
    double mv = molarVolume();
    double sqt = sqrt(temperature());
    double vpb = mv + m_b_current;
    double vmb = mv - m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    double pres = pressure();

    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = (- RT() * log(pres * mv / RT())
                 + RT() * log(mv / vmb)
                 + RT() * b_vec_Curr_[k] / vmb
                 - 2.0 * m_pp[k] / (m_b_current * sqt) * log(vpb/mv)
                 + m_a_current * b_vec_Curr_[k] / (m_b_current * m_b_current * sqt) * log(vpb/mv)
                 - m_a_current / (m_b_current * sqt) * (b_vec_Curr_[k]/vpb)
                );
    }
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(ac[k]/RT());
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void RedlichKwongMFTP::getChemPotentials(span<double> mu) const
{
    getGibbs_ref(mu);
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT()*(log(xx));
    }

    double mv = molarVolume();
    double sqt = sqrt(temperature());
    double vpb = mv + m_b_current;
    double vmb = mv - m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    double pres = pressure();
    double refP = refPressure();

    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += (RT() * log(pres/refP) - RT() * log(pres * mv / RT())
                  + RT() * log(mv / vmb)
                  + RT() * b_vec_Curr_[k] / vmb
                  - 2.0 * m_pp[k] / (m_b_current * sqt) * log(vpb/mv)
                  + m_a_current * b_vec_Curr_[k] / (m_b_current * m_b_current * sqt) * log(vpb/mv)
                  - m_a_current / (m_b_current * sqt) * (b_vec_Curr_[k]/vpb)
                 );
    }
}

void RedlichKwongMFTP::getPartialMolarEnthalpies(span<double> hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar.begin(), hbar.end(), hbar.begin(), RT());

    // We calculate dpdni_
    double TKelvin = temperature();
    double mv = molarVolume();
    double sqt = sqrt(TKelvin);
    double vpb = mv + m_b_current;
    double vmb = mv - m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    for (size_t k = 0; k < m_kk; k++) {
        dpdni_[k] = RT()/vmb + RT() * b_vec_Curr_[k] / (vmb * vmb) - 2.0 * m_pp[k] / (sqt * mv * vpb)
                    + m_a_current * b_vec_Curr_[k]/(sqt * mv * vpb * vpb);
    }
    double dadt = da_dt();
    double fac = TKelvin * dadt - 3.0 * m_a_current / 2.0;

    for (size_t k = 0; k < m_kk; k++) {
        m_workS[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_workS[k] += 2.0 * moleFractions_[i] * TKelvin * a_coeff_vec(1,counter) - 3.0 * moleFractions_[i] * a_vec_Curr_[counter];
        }
    }

    pressureDerivatives();
    double fac2 = mv + TKelvin * dpdT_ / dpdV_;
    for (size_t k = 0; k < m_kk; k++) {
        double hE_v = (mv * dpdni_[k] - RT() - b_vec_Curr_[k]/ (m_b_current * m_b_current * sqt) * log(vpb/mv)*fac
                       + 1.0 / (m_b_current * sqt) * log(vpb/mv) * m_workS[k]
                       +  b_vec_Curr_[k] / vpb / (m_b_current * sqt) * fac);
        hbar[k] = hbar[k] + hE_v;
        hbar[k] -= fac2 * dpdni_[k];
    }
}

void RedlichKwongMFTP::getPartialMolarEntropies(span<double> sbar) const
{
    getEntropy_R_ref(sbar);
    scale(sbar.begin(), sbar.end(), sbar.begin(), GasConstant);
    double TKelvin = temperature();
    double sqt = sqrt(TKelvin);
    double mv = molarVolume();
    double refP = refPressure();

    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (- log(xx));
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_workS[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_workS[k] += moleFractions_[i] * a_coeff_vec(1,counter);
        }
    }

    double dadt = da_dt();
    double fac = dadt - m_a_current / (2.0 * TKelvin);
    double vmb = mv - m_b_current;
    double vpb = mv + m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -=(GasConstant * log(GasConstant * TKelvin / (refP * mv))
                   + GasConstant
                   + GasConstant * log(mv/vmb)
                   + GasConstant * b_vec_Curr_[k]/vmb
                   + m_pp[k]/(m_b_current * TKelvin * sqt) * log(vpb/mv)
                   - 2.0 * m_workS[k]/(m_b_current * sqt) * log(vpb/mv)
                   + b_vec_Curr_[k] / (m_b_current * m_b_current * sqt) * log(vpb/mv) * fac
                   - 1.0 / (m_b_current * sqt) * b_vec_Curr_[k] / vpb * fac
                  );
    }

    pressureDerivatives();
    getPartialMolarVolumes(m_partialMolarVolumes);
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= -m_partialMolarVolumes[k] * dpdT_;
    }
}

void RedlichKwongMFTP::getPartialMolarIntEnergies(span<double> ubar) const
{
    // u_k = h_k - P * v_k
    getPartialMolarVolumes(m_partialMolarVolumes);
    getPartialMolarEnthalpies(ubar);
    double p = pressure();
    for (size_t k = 0; k < nSpecies(); k++) {
        ubar[k] -= p * m_partialMolarVolumes[k];
    }
}

void RedlichKwongMFTP::getPartialMolarVolumes(span<double> vbar) const
{
    checkArraySize("RedlichKwongMFTP::getPartialMolarVolumes", vbar.size(), m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_workS[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_workS[k] += moleFractions_[i] * a_coeff_vec(1,counter);
        }
    }

    double sqt = sqrt(temperature());
    double mv = molarVolume();
    double vmb = mv - m_b_current;
    double vpb = mv + m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
        double num = (RT() + RT() * m_b_current/ vmb + RT() * b_vec_Curr_[k] / vmb
                          + RT() * m_b_current * b_vec_Curr_[k] /(vmb * vmb)
                          - 2.0 * m_pp[k] / (sqt * vpb)
                          + m_a_current * b_vec_Curr_[k] / (sqt * vpb * vpb)
                         );
        double denom = (pressure() + RT() * m_b_current/(vmb * vmb) - m_a_current / (sqt * vpb * vpb)
                           );
        vbar[k] = num / denom;
    }
}

bool RedlichKwongMFTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        a_vec_Curr_.resize(m_kk * m_kk, 0.0);

        // Initialize a_vec and b_vec to NaN, to screen for species with
        //     pureFluidParameters which are undefined in the input file:
        b_vec_Curr_.push_back(NAN);
        a_coeff_vec.resize(2, m_kk * m_kk, NAN);

        m_pp.push_back(0.0);
        m_coeffSource.push_back(CoeffSource::EoS);
        m_partialMolarVolumes.push_back(0.0);
        dpdni_.push_back(0.0);
    }
    return added;
}

void RedlichKwongMFTP::initThermo()
{
    // Contents of 'critical-properties.yaml', loaded later if needed
    AnyMap critPropsDb;
    std::unordered_map<string, AnyMap*> dbSpecies;

    for (auto& [name, species] : m_species) {
        auto& data = species->input;
        size_t k = speciesIndex(name, true);
        if (!isnan(a_coeff_vec(0, k + m_kk * k))) {
            continue;
        }
        bool foundCoeffs = false;

        if (data.hasKey("equation-of-state") &&
            data["equation-of-state"].hasMapWhere("model", "Redlich-Kwong"))
        {
            // Read a and b coefficients from species 'input' information (that is, as
            // specified in a YAML input file) specific to the Redlich-Kwong model
            auto eos = data["equation-of-state"].getMapWhere(
                "model", "Redlich-Kwong");

            if (eos.hasKey("a") && eos.hasKey("b")) {
                double a0 = 0, a1 = 0;
                if (eos["a"].isScalar()) {
                    a0 = eos.convert("a", "Pa*m^6/kmol^2*K^0.5");
                } else {
                    auto avec = eos["a"].asVector<AnyValue>(2);
                    a0 = eos.units().convert(avec[0], "Pa*m^6/kmol^2*K^0.5");
                    a1 = eos.units().convert(avec[1], "Pa*m^6/kmol^2/K^0.5");
                }
                double b = eos.convert("b", "m^3/kmol");
                foundCoeffs = true;
                setSpeciesCoeffs(name, a0, a1, b);
                m_coeffSource[k] = CoeffSource::EoS;
            }

            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& [name2, coeff] : binary_a) {
                    double a0 = 0, a1 = 0;
                    if (coeff.isScalar()) {
                        a0 = units.convert(coeff, "Pa*m^6/kmol^2*K^0.5");
                    } else {
                        auto avec = coeff.asVector<AnyValue>(2);
                        a0 = units.convert(avec[0], "Pa*m^6/kmol^2*K^0.5");
                        a1 = units.convert(avec[1], "Pa*m^6/kmol^2/K^0.5");
                    }
                    setBinaryCoeffs(name, name2, a0, a1);
                }
            }

            if (foundCoeffs) {
                continue;
            }
        }

        // Coefficients have not been populated from model-specific input
        double Tc = NAN, Pc = NAN;
        if (data.hasKey("critical-parameters")) {
            // Use critical state information stored in the species entry to
            // calculate a and b
            auto& critProps = data["critical-parameters"].as<AnyMap>();
            Tc = critProps.convert("critical-temperature", "K");
            Pc = critProps.convert("critical-pressure", "Pa");
            m_coeffSource[k] = CoeffSource::CritProps;
        } else {
            // Search 'crit-properties.yaml' to find Tc and Pc. Load data if needed
            if (critPropsDb.empty()) {
                critPropsDb = AnyMap::fromYamlFile("critical-properties.yaml");
                dbSpecies = critPropsDb["species"].asMap("name");
            }

            // All names in critical-properties.yaml are upper case
            auto ucName = boost::algorithm::to_upper_copy(name);
            if (dbSpecies.count(ucName)) {
                auto& spec = *dbSpecies.at(ucName);
                auto& critProps = spec["critical-parameters"].as<AnyMap>();
                Tc = critProps.convert("critical-temperature", "K");
                Pc = critProps.convert("critical-pressure", "Pa");
                m_coeffSource[k] = CoeffSource::Database;
            }
        }

        // Check if critical properties were found in either location
        if (!isnan(Tc)) {
            // Assuming no temperature dependence (that is, a1 = 0)
            double a = omega_a * pow(GasConstant, 2) * pow(Tc, 2.5) / Pc;
            double b = omega_b * GasConstant * Tc / Pc;
            setSpeciesCoeffs(name, a, 0.0, b);
        } else {
            throw InputFileError("RedlichKwongMFTP::initThermo", data,
                    "No critical property or Redlich-Kwong parameters found "
                    "for species {}.", name);
        }
    }
}

void RedlichKwongMFTP::getSpeciesParameters(const string& name,
                                            AnyMap& speciesNode) const
{
    MixtureFugacityTP::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name, true);
    if (m_coeffSource[k] == CoeffSource::EoS) {
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Redlich-Kwong", true);

        size_t counter = k + m_kk * k;
        if (a_coeff_vec(1, counter) != 0.0) {
            vector<AnyValue> coeffs(2);
            coeffs[0].setQuantity(a_coeff_vec(0, counter), "Pa*m^6/kmol^2*K^0.5");
            coeffs[1].setQuantity(a_coeff_vec(1, counter), "Pa*m^6/kmol^2/K^0.5");
            eosNode["a"] = std::move(coeffs);
        } else {
            eosNode["a"].setQuantity(a_coeff_vec(0, counter),
                                    "Pa*m^6/kmol^2*K^0.5");
        }
        eosNode["b"].setQuantity(b_vec_Curr_[k], "m^3/kmol");
    } else if (m_coeffSource[k] == CoeffSource::CritProps) {
        auto& critProps = speciesNode["critical-parameters"];
        double a = a_coeff_vec(0, k + m_kk * k);
        double b = b_vec_Curr_[k];
        double Tc = pow(a * omega_b / (b * omega_a * GasConstant), 2.0/3.0);
        double Pc = omega_b * GasConstant * Tc / b;
        critProps["critical-temperature"].setQuantity(Tc, "K");
        critProps["critical-pressure"].setQuantity(Pc, "Pa");
    }

    if (m_binaryParameters.count(name)) {
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Redlich-Kwong", true);
        AnyMap bin_a;
        for (const auto& [name2, coeffs] : m_binaryParameters.at(name)) {
            if (coeffs.second == 0) {
                bin_a[name2].setQuantity(coeffs.first, "Pa*m^6/kmol^2*K^0.5");
            } else {
                vector<AnyValue> C(2);
                C[0].setQuantity(coeffs.first, "Pa*m^6/kmol^2*K^0.5");
                C[1].setQuantity(coeffs.second, "Pa*m^6/kmol^2/K^0.5");
                bin_a[name2] = std::move(C);
            }
        }
        eosNode["binary-a"] = std::move(bin_a);
    }
}

double RedlichKwongMFTP::sresid() const
{
    // note this agrees with tpx
    double rho = density();
    double mmw = meanMolecularWeight();
    double molarV = mmw / rho;
    double hh = m_b_current / molarV;
    double zz = z();
    double dadt = da_dt();
    double T = temperature();
    double sqT = sqrt(T);
    double fac = dadt - m_a_current / (2.0 * T);
    double sresid_mol_R = log(zz*(1.0 - hh)) + log(1.0 + hh) * fac / (sqT * GasConstant * m_b_current);
    return GasConstant * sresid_mol_R;
}

double RedlichKwongMFTP::hresid() const
{
    // note this agrees with tpx
    double rho = density();
    double mmw = meanMolecularWeight();
    double molarV = mmw / rho;
    double hh = m_b_current / molarV;
    double zz = z();
    double dadt = da_dt();
    double T = temperature();
    double sqT = sqrt(T);
    double fac = T * dadt - 3.0 *m_a_current / (2.0);
    return GasConstant * T * (zz - 1.0) + fac * log(1.0 + hh) / (sqT * m_b_current);
}

double RedlichKwongMFTP::liquidVolEst(double TKelvin, double& presGuess) const
{
    double v = m_b_current * 1.1;
    double atmp;
    double btmp;
    calculateAB(TKelvin, atmp, btmp);
    double pres = std::max(psatEst(TKelvin), presGuess);
    double Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
        int nsol = solveCubic(TKelvin, pres, atmp, btmp, Vroot);
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

double RedlichKwongMFTP::densityCalc(double TKelvin, double presPa, int phaseRequested,
                                     double rhoguess)
{
    // It's necessary to set the temperature so that m_a_current is set correctly.
    setTemperature(TKelvin);
    double tcrit = critTemperature();
    double mmw = meanMolecularWeight();
    if (rhoguess == -1.0) {
        if (phaseRequested >= FLUID_LIQUID_0) {
                    double lqvol = liquidVolEst(TKelvin, presPa);
                    rhoguess = mmw / lqvol;
                }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to the routine
            rhoguess = presPa * mmw / (GasConstant * TKelvin);
    }


    double volguess = mmw / rhoguess;
    NSolns_ = solveCubic(TKelvin, presPa, m_a_current, m_b_current, Vroot_);

    // Ensure root ordering used for branch selection only contains physical roots.
    double vmin = std::max(0.0, m_b_current * (1.0 + 1e-12));
    vector<double> physicalRoots;
    for (double root : Vroot_) {
        if (std::isfinite(root) && root > vmin) {
            physicalRoots.push_back(root);
        }
    }
    std::sort(physicalRoots.begin(), physicalRoots.end());
    if (physicalRoots.empty()) {
        return -1.0;
    } else if (physicalRoots.size() == 1) {
        // Preserve branch-request semantics for single-root states:
        // return -2 when the requested phase is inconsistent with the
        // single branch identified by solveCubic() sign convention.
        if ((phaseRequested == FLUID_GAS && NSolns_ < 0)
            || (phaseRequested >= FLUID_LIQUID_0 && NSolns_ > 0))
        {
            return -2.0;
        }
        // Otherwise, accept the physically admissible root directly.
        return mmw / physicalRoots[0];
    } else if (physicalRoots.size() == 2) {
        Vroot_[0] = physicalRoots[0];
        Vroot_[1] = 0.5 * (physicalRoots[0] + physicalRoots[1]);
        Vroot_[2] = physicalRoots[1];
    } else {
        Vroot_[0] = physicalRoots[0];
        Vroot_[1] = physicalRoots[1];
        Vroot_[2] = physicalRoots[2];
    }

    double molarVolLast = Vroot_[0];
    if (NSolns_ >= 2) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            molarVolLast = Vroot_[0];
        } else if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[2];
        } else {
            if (volguess > Vroot_[1]) {
                molarVolLast = Vroot_[2];
            } else {
                molarVolLast = Vroot_[0];
            }
        }
    } else if (NSolns_ == 1) {
        if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT || phaseRequested == FLUID_UNDEFINED) {
            molarVolLast = Vroot_[0];
        } else {
            return -2.0;
        }
    } else if (NSolns_ == -1) {
        if (phaseRequested >= FLUID_LIQUID_0 || phaseRequested == FLUID_UNDEFINED || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[0];
        } else if (TKelvin > tcrit) {
            molarVolLast = Vroot_[0];
        } else {
            return -2.0;
        }
    } else {
        molarVolLast = Vroot_[0];
        return -1.0;
    }
    return mmw / molarVolLast;
}

double RedlichKwongMFTP::dpdVCalc(double TKelvin, double molarVol, double& presCalc) const
{
    double sqt = sqrt(TKelvin);
    presCalc = GasConstant * TKelvin / (molarVol - m_b_current)
               - m_a_current / (sqt * molarVol * (molarVol + m_b_current));

    double vpb = molarVol + m_b_current;
    double vmb = molarVol - m_b_current;
    double dpdv = (- GasConstant * TKelvin / (vmb * vmb)
                       + m_a_current * (2 * molarVol + m_b_current) / (sqt * molarVol * molarVol * vpb * vpb));
    return dpdv;
}

double RedlichKwongMFTP::isothermalCompressibility() const
{
    pressureDerivatives();
    return -1 / (molarVolume() * dpdV_);
}

double RedlichKwongMFTP::thermalExpansionCoeff() const
{
    pressureDerivatives();
    return -dpdT_ / (molarVolume() * dpdV_);
}

double RedlichKwongMFTP::soundSpeed() const
{
    pressureDerivatives();
    return molarVolume() * sqrt(-cp_mole() / cv_mole() * dpdV_ / meanMolecularWeight());
}

void RedlichKwongMFTP::pressureDerivatives() const
{
    double TKelvin = temperature();
    double mv = molarVolume();
    double pres;

    dpdV_ = dpdVCalc(TKelvin, mv, pres);
    double sqt = sqrt(TKelvin);
    double vpb = mv + m_b_current;
    double vmb = mv - m_b_current;
    double dadt = da_dt();
    double fac = dadt - m_a_current/(2.0 * TKelvin);
    dpdT_ = (GasConstant / vmb - fac / (sqt * mv * vpb));
}

void RedlichKwongMFTP::updateMixingExpressions()
{
    double temp = temperature();
    if (m_formTempParam == 1) {
        for (size_t i = 0; i < m_kk; i++) {
            for (size_t j = 0; j < m_kk; j++) {
                size_t counter = i * m_kk + j;
                a_vec_Curr_[counter] = a_coeff_vec(0,counter) + a_coeff_vec(1,counter) * temp;
            }
        }
    }

    m_b_current = 0.0;
    m_a_current = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        m_b_current += moleFractions_[i] * b_vec_Curr_[i];
        for (size_t j = 0; j < m_kk; j++) {
            m_a_current += a_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
        }
    }
    if (isnan(m_b_current)) {
        // One or more species do not have specified coefficients.
        fmt::memory_buffer b;
        for (size_t k = 0; k < m_kk; k++) {
            if (isnan(b_vec_Curr_[k])) {
                if (b.size() > 0) {
                    fmt_append(b, ", {}", speciesName(k));
                } else {
                    fmt_append(b, "{}", speciesName(k));
                }
            }
        }
        throw CanteraError("RedlichKwongMFTP::updateMixingExpressions",
            "Missing Redlich-Kwong coefficients for species: {}", to_string(b));
    }
}

void RedlichKwongMFTP::calculateAB(double temp, double& aCalc, double& bCalc) const
{
    bCalc = 0.0;
    aCalc = 0.0;
    if (m_formTempParam == 1) {
        for (size_t i = 0; i < m_kk; i++) {
            bCalc += moleFractions_[i] * b_vec_Curr_[i];
            for (size_t j = 0; j < m_kk; j++) {
                size_t counter = i * m_kk + j;
                double a_vec_Curr = a_coeff_vec(0,counter) + a_coeff_vec(1,counter) * temp;
                aCalc += a_vec_Curr * moleFractions_[i] * moleFractions_[j];
            }
        }
    } else {
        for (size_t i = 0; i < m_kk; i++) {
            bCalc += moleFractions_[i] * b_vec_Curr_[i];
            for (size_t j = 0; j < m_kk; j++) {
                size_t counter = i * m_kk + j;
                double a_vec_Curr = a_coeff_vec(0,counter);
                aCalc += a_vec_Curr * moleFractions_[i] * moleFractions_[j];
            }
        }
    }
}

double RedlichKwongMFTP::da_dt() const
{
    double dadT = 0.0;
    if (m_formTempParam == 1) {
        for (size_t i = 0; i < m_kk; i++) {
            for (size_t j = 0; j < m_kk; j++) {
                size_t counter = i * m_kk + j;
                dadT+= a_coeff_vec(1,counter) * moleFractions_[i] * moleFractions_[j];
            }
        }
    }
    return dadT;
}

void RedlichKwongMFTP::calcCriticalConditions(double& pc, double& tc, double& vc) const
{
    double a0 = 0.0;
    double aT = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j <m_kk; j++) {
            size_t counter = i + m_kk * j;
            a0 += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(0, counter);
            aT += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(1, counter);
        }
    }
    double a = m_a_current;
    double b = m_b_current;
    if (m_formTempParam != 0) {
        a = a0;
    }
    if (b <= 0.0) {
        tc = 1000000.;
        pc = 1.0E13;
        vc = omega_vc * GasConstant * tc / pc;
        return;
    }
    if (a <= 0.0) {
        tc = 0.0;
        pc = 0.0;
        vc = 2.0 * b;
        return;
    }
    double tmp = a * omega_b / (b * omega_a * GasConstant);
    double pp = 2./3.;
    double sqrttc, f, dfdt, deltatc;

    if (m_formTempParam == 0) {
        tc = pow(tmp, pp);
    } else {
        tc = pow(tmp, pp);
        for (int j = 0; j < 10; j++) {
            sqrttc = sqrt(tc);
            f = omega_a * b * GasConstant * tc * sqrttc / omega_b - aT * tc - a0;
            dfdt = 1.5 * omega_a * b * GasConstant * sqrttc / omega_b - aT;
            deltatc = - f / dfdt;
            tc += deltatc;
        }
        if (deltatc > 0.1) {
            throw CanteraError("RedlichKwongMFTP::calcCriticalConditions",
                "didn't converge");
        }
    }

    pc = omega_b * GasConstant * tc / b;
    vc = omega_vc * GasConstant * tc / pc;
}

int RedlichKwongMFTP::solveCubic(double T, double pres, double a, double b,
                                 span<double> Vroot) const
{

    // Derive the coefficients of the cubic polynomial to solve.
    double an = 1.0;
    double bn = - GasConstant * T / pres;
    double sqt = sqrt(T);
    double cn = - (GasConstant * T * b / pres - a/(pres * sqt) + b * b);
    double dn = - (a * b / (pres * sqt));

    double tmp = a * omega_b / (b * omega_a * GasConstant);
    double pp = 2./3.;
    double tc = pow(tmp, pp);
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    int nSolnValues = MixtureFugacityTP::solveCubic(T, pres, a, b, a, Vroot, an, bn, cn, dn, tc, vc);

    return nSolnValues;
}

}
