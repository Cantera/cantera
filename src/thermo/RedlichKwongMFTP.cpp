//! @file RedlichKwongMFTP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <boost/math/tools/roots.hpp>

#include <algorithm>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

const doublereal RedlichKwongMFTP::omega_a = 4.27480233540E-01;
const doublereal RedlichKwongMFTP::omega_b = 8.66403499650E-02;
const doublereal RedlichKwongMFTP::omega_vc = 3.33333333333333E-01;

RedlichKwongMFTP::RedlichKwongMFTP(const std::string& infile, const std::string& id_) :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    fill_n(Vroot_, 3, 0.0);
    initThermoFile(infile, id_);
}

RedlichKwongMFTP::RedlichKwongMFTP(XML_Node& phaseRefRoot, const std::string& id_) :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    fill_n(Vroot_, 3, 0.0);
    importPhase(phaseRefRoot, this);
}

void RedlichKwongMFTP::setSpeciesCoeffs(const std::string& species,
                                        double a0, double a1, double b)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
        throw CanteraError("RedlichKwongMFTP::setSpeciesCoeffs",
            "Unknown species '{}'.", species);
    }

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
    a_coeff_vec.getRow(0, a_vec_Curr_.data());
    b_vec_Curr_[k] = b;
}

void RedlichKwongMFTP::setBinaryCoeffs(const std::string& species_i,
        const std::string& species_j, double a0, double a1)
{
    size_t ki = speciesIndex(species_i);
    if (ki == npos) {
        throw CanteraError("RedlichKwongMFTP::setBinaryCoeffs",
            "Unknown species '{}'.", species_i);
    }
    size_t kj = speciesIndex(species_j);
    if (kj == npos) {
        throw CanteraError("RedlichKwongMFTP::setBinaryCoeffs",
            "Unknown species '{}'.", species_j);
    }

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

doublereal RedlichKwongMFTP::cp_mole() const
{
    _updateReferenceStateThermo();
    doublereal TKelvin = temperature();
    doublereal sqt = sqrt(TKelvin);
    doublereal mv = molarVolume();
    doublereal vpb = mv + m_b_current;
    pressureDerivatives();
    doublereal cpref = GasConstant * mean_X(m_cp0_R);
    doublereal dadt = da_dt();
    doublereal fac = TKelvin * dadt - 3.0 * m_a_current / 2.0;
    doublereal dHdT_V = (cpref + mv * dpdT_ - GasConstant - 1.0 / (2.0 * m_b_current * TKelvin * sqt) * log(vpb/mv) * fac
                         +1.0/(m_b_current * sqt) * log(vpb/mv) * (-0.5 * dadt));
    return dHdT_V - (mv + TKelvin * dpdT_ / dpdV_) * dpdT_;
}

doublereal RedlichKwongMFTP::cv_mole() const
{
    _updateReferenceStateThermo();
    doublereal TKelvin = temperature();
    doublereal sqt = sqrt(TKelvin);
    doublereal mv = molarVolume();
    doublereal vpb = mv + m_b_current;
    doublereal cvref = GasConstant * (mean_X(m_cp0_R) - 1.0);
    doublereal dadt = da_dt();
    doublereal fac = TKelvin * dadt - 3.0 * m_a_current / 2.0;
    return (cvref - 1.0/(2.0 * m_b_current * TKelvin * sqt) * log(vpb/mv)*fac
            +1.0/(m_b_current * sqt) * log(vpb/mv)*(-0.5*dadt));
}

doublereal RedlichKwongMFTP::pressure() const
{
    _updateReferenceStateThermo();

    //  Get a copy of the private variables stored in the State object
    doublereal T = temperature();
    double molarV = meanMolecularWeight() / density();
    double pp = GasConstant * T/(molarV - m_b_current) - m_a_current/(sqrt(T) * molarV * (molarV + m_b_current));
    return pp;
}

doublereal RedlichKwongMFTP::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void RedlichKwongMFTP::getActivityCoefficients(doublereal* ac) const
{
    doublereal mv = molarVolume();
    doublereal sqt = sqrt(temperature());
    doublereal vpb = mv + m_b_current;
    doublereal vmb = mv - m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    doublereal pres = pressure();

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

void RedlichKwongMFTP::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= 1.0 / RT();
    }
}

void RedlichKwongMFTP::getChemPotentials(doublereal* mu) const
{
    getGibbs_ref(mu);
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT()*(log(xx));
    }

    doublereal mv = molarVolume();
    doublereal sqt = sqrt(temperature());
    doublereal vpb = mv + m_b_current;
    doublereal vmb = mv - m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    doublereal pres = pressure();
    doublereal refP = refPressure();

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

void RedlichKwongMFTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar, hbar+m_kk, hbar, RT());

    // We calculate dpdni_
    doublereal TKelvin = temperature();
    doublereal mv = molarVolume();
    doublereal sqt = sqrt(TKelvin);
    doublereal vpb = mv + m_b_current;
    doublereal vmb = mv - m_b_current;
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
    doublereal dadt = da_dt();
    doublereal fac = TKelvin * dadt - 3.0 * m_a_current / 2.0;

    for (size_t k = 0; k < m_kk; k++) {
        m_tmpV[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_tmpV[k] += 2.0 * moleFractions_[i] * TKelvin * a_coeff_vec(1,counter) - 3.0 * moleFractions_[i] * a_vec_Curr_[counter];
        }
    }

    pressureDerivatives();
    doublereal fac2 = mv + TKelvin * dpdT_ / dpdV_;
    for (size_t k = 0; k < m_kk; k++) {
        double hE_v = (mv * dpdni_[k] - RT() - b_vec_Curr_[k]/ (m_b_current * m_b_current * sqt) * log(vpb/mv)*fac
                       + 1.0 / (m_b_current * sqt) * log(vpb/mv) * m_tmpV[k]
                       +  b_vec_Curr_[k] / vpb / (m_b_current * sqt) * fac);
        hbar[k] = hbar[k] + hE_v;
        hbar[k] -= fac2 * dpdni_[k];
    }
}

void RedlichKwongMFTP::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R_ref(sbar);
    scale(sbar, sbar+m_kk, sbar, GasConstant);
    doublereal TKelvin = temperature();
    doublereal sqt = sqrt(TKelvin);
    doublereal mv = molarVolume();
    doublereal refP = refPressure();

    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
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
        m_tmpV[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_tmpV[k] += moleFractions_[i] * a_coeff_vec(1,counter);
        }
    }

    doublereal dadt = da_dt();
    doublereal fac = dadt - m_a_current / (2.0 * TKelvin);
    doublereal vmb = mv - m_b_current;
    doublereal vpb = mv + m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -=(GasConstant * log(GasConstant * TKelvin / (refP * mv))
                   + GasConstant
                   + GasConstant * log(mv/vmb)
                   + GasConstant * b_vec_Curr_[k]/vmb
                   + m_pp[k]/(m_b_current * TKelvin * sqt) * log(vpb/mv)
                   - 2.0 * m_tmpV[k]/(m_b_current * sqt) * log(vpb/mv)
                   + b_vec_Curr_[k] / (m_b_current * m_b_current * sqt) * log(vpb/mv) * fac
                   - 1.0 / (m_b_current * sqt) * b_vec_Curr_[k] / vpb * fac
                  );
    }

    pressureDerivatives();
    getPartialMolarVolumes(m_partialMolarVolumes.data());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= -m_partialMolarVolumes[k] * dpdT_;
    }
}

void RedlichKwongMFTP::getPartialMolarIntEnergies(doublereal* ubar) const
{
    // u_k = h_k - P * v_k
    getPartialMolarVolumes(m_partialMolarVolumes.data());
    getPartialMolarEnthalpies(ubar);
    double p = pressure();
    for (size_t k = 0; k < nSpecies(); k++) {
        ubar[k] -= p * m_partialMolarVolumes[k];
    }
}

void RedlichKwongMFTP::getPartialMolarVolumes(doublereal* vbar) const
{
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * a_vec_Curr_[counter];
        }
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_tmpV[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_tmpV[k] += moleFractions_[i] * a_coeff_vec(1,counter);
        }
    }

    doublereal sqt = sqrt(temperature());
    doublereal mv = molarVolume();
    doublereal vmb = mv - m_b_current;
    doublereal vpb = mv + m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
        doublereal num = (RT() + RT() * m_b_current/ vmb + RT() * b_vec_Curr_[k] / vmb
                          + RT() * m_b_current * b_vec_Curr_[k] /(vmb * vmb)
                          - 2.0 * m_pp[k] / (sqt * vpb)
                          + m_a_current * b_vec_Curr_[k] / (sqt * vpb * vpb)
                         );
        doublereal denom = (pressure() + RT() * m_b_current/(vmb * vmb) - m_a_current / (sqt * vpb * vpb)
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

void RedlichKwongMFTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thermoNode = phaseNode.child("thermo");
        std::string model = thermoNode["model"];
        if (model != "RedlichKwong" && model != "RedlichKwongMFTP") {
            throw CanteraError("RedlichKwongMFTP::initThermoXML",
                               "Unknown thermo model : " + model);
        }

        // Reset any coefficients which may have been set using values from
        // 'critical-properties.yaml' as part of non-XML initialization, so that
        // off-diagonal elements can be correctly initialized
        a_coeff_vec.data().assign(a_coeff_vec.data().size(), NAN);

        // Go get all of the coefficients and factors in the
        // activityCoefficients XML block
        if (thermoNode.hasChild("activityCoefficients")) {
            XML_Node& acNode = thermoNode.child("activityCoefficients");

            // Loop through the children and read out fluid parameters.  Process
            //   all the pureFluidParameters, first:
            // Loop back through the "activityCoefficients" children and process the
            // crossFluidParameters in the XML tree:
            for (size_t i = 0; i < acNode.nChildren(); i++) {
                XML_Node& xmlACChild = acNode.child(i);
                if (caseInsensitiveEquals(xmlACChild.name(), "purefluidparameters")) {
                    readXMLPureFluid(xmlACChild);
                } else if (caseInsensitiveEquals(xmlACChild.name(), "crossfluidparameters")) {
                    readXMLCrossFluid(xmlACChild);
                }
            }
        }
        // If any species exist which have undefined pureFluidParameters,
        // search the database in 'critical-properties.yaml' to find critical
        // temperature and pressure to calculate a and b.

        // Loop through all species in the CTI file
        size_t iSpecies = 0;

        for (size_t i = 0; i < m_kk; i++) {
            string iName = speciesName(i);

            // Get the index of the species
            iSpecies = speciesIndex(iName);

            // Check if a and b are already populated (only the diagonal elements of a).
            size_t counter = iSpecies + m_kk * iSpecies;

            // If not, then search the database:
            if (isnan(a_coeff_vec(0, counter))) {

                vector<double> coeffArray;

                // Search the database for the species name and calculate
                // coefficients a and b, from critical properties:
                // coeffArray[0] = a0, coeffArray[1] = b;
                coeffArray = getCoeff(iName);

                // Check if species was found in the database of critical properties,
                // and assign the results
                if (!isnan(coeffArray[0])) {
                    //Assuming no temperature dependence (i,e a1 = 0)
                    setSpeciesCoeffs(iName, coeffArray[0], 0.0, coeffArray[1]);
                    m_coeffSource[i] = CoeffSource::EoS;
                }
            }
        }
    }

    MixtureFugacityTP::initThermoXML(phaseNode, id);
}

void RedlichKwongMFTP::initThermo()
{
    // Contents of 'critical-properties.yaml', loaded later if needed
    AnyMap critPropsDb;
    std::unordered_map<std::string, AnyMap*> dbSpecies;

    for (auto& item : m_species) {
        auto& data = item.second->input;
        size_t k = speciesIndex(item.first);
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
                setSpeciesCoeffs(item.first, a0, a1, b);
                m_coeffSource[k] = CoeffSource::EoS;
            }

            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& item2 : binary_a) {
                    double a0 = 0, a1 = 0;
                    if (item2.second.isScalar()) {
                        a0 = units.convert(item2.second, "Pa*m^6/kmol^2*K^0.5");
                    } else {
                        auto avec = item2.second.asVector<AnyValue>(2);
                        a0 = units.convert(avec[0], "Pa*m^6/kmol^2*K^0.5");
                        a1 = units.convert(avec[1], "Pa*m^6/kmol^2/K^0.5");
                    }
                    setBinaryCoeffs(item.first, item2.first, a0, a1);
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
            auto ucName = boost::algorithm::to_upper_copy(item.first);
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
            setSpeciesCoeffs(item.first, a, 0.0, b);
        }
        // @todo: After XML is removed, this can throw an InputFileError if neither
        // R-K parameters or critical properties were found.
    }
}

void RedlichKwongMFTP::getSpeciesParameters(const std::string& name,
                                            AnyMap& speciesNode) const
{
    MixtureFugacityTP::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name);
    checkSpeciesIndex(k);
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
        for (const auto& item : m_binaryParameters.at(name)) {
            if (item.second.second == 0) {
                bin_a[item.first].setQuantity(item.second.first, "Pa*m^6/kmol^2*K^0.5");
            } else {
                vector<AnyValue> coeffs(2);
                coeffs[0].setQuantity(item.second.first, "Pa*m^6/kmol^2*K^0.5");
                coeffs[1].setQuantity(item.second.second, "Pa*m^6/kmol^2/K^0.5");
                bin_a[item.first] = std::move(coeffs);
            }
        }
        eosNode["binary-a"] = std::move(bin_a);
    }
}

vector<double> RedlichKwongMFTP::getCoeff(const std::string& iName)
{
    warn_deprecated("RedlichKwongMFTP::getCoeff", "To be removed after Cantera 2.6. "
                    "Use of critical-properties.yaml is integrated into initThermo() "
                    "for YAML input files.");
    vector_fp spCoeff{NAN, NAN};
    AnyMap data = AnyMap::fromYamlFile("critical-properties.yaml");
    const auto& species = data["species"].asMap("name");

    // All names in critical-properties.yaml are upper case
    auto ucName = boost::algorithm::to_upper_copy(iName);
    if (species.count(ucName)) {
        auto& critProps = species.at(ucName)->at("critical-parameters").as<AnyMap>();
        double Tc = critProps.convert("critical-temperature", "K");
        double Pc = critProps.convert("critical-pressure", "Pa");

        // calculate pure fluid parameters a_k and b_k based on T_c and P_c assuming
        // no temperature dependence
        spCoeff[0] = omega_a * pow(GasConstant, 2) * pow(Tc, 2.5) / Pc; //coeff a
        spCoeff[1] = omega_b * GasConstant * Tc / Pc; // coeff b
    }
    return spCoeff;
}

void RedlichKwongMFTP::readXMLPureFluid(XML_Node& pureFluidParam)
{
    string xname = pureFluidParam.name();
    if (xname != "pureFluidParameters") {
        throw CanteraError("RedlichKwongMFTP::readXMLPureFluid",
                           "Incorrect name for processing this routine: " + xname);
    }

    double a0 = 0.0;
    double a1 = 0.0;
    double b = 0.0;
    for (size_t iChild = 0; iChild < pureFluidParam.nChildren(); iChild++) {
        XML_Node& xmlChild = pureFluidParam.child(iChild);
        string nodeName = toLowerCopy(xmlChild.name());

        if (nodeName == "a_coeff") {
            vector_fp vParams;
            string iModel = toLowerCopy(xmlChild.attrib("model"));
            getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");

            if (iModel == "constant" && vParams.size() == 1) {
                a0 = vParams[0];
                a1 = 0;
            } else if (iModel == "linear_a" && vParams.size() == 2) {
                a0 = vParams[0];
                a1 = vParams[1];
            } else {
                throw CanteraError("RedlichKwongMFTP::readXMLPureFluid",
                    "unknown model or incorrect number of parameters");
            }

        } else if (nodeName == "b_coeff") {
            b = getFloatCurrent(xmlChild, "toSI");
        }
    }
    setSpeciesCoeffs(pureFluidParam.attrib("species"), a0, a1, b);
    m_coeffSource[speciesIndex(pureFluidParam.attrib("species"))] = CoeffSource::EoS;
}

void RedlichKwongMFTP::readXMLCrossFluid(XML_Node& CrossFluidParam)
{
    string xname = CrossFluidParam.name();
    if (xname != "crossFluidParameters") {
        throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid",
                           "Incorrect name for processing this routine: " + xname);
    }

    string iName = CrossFluidParam.attrib("species1");
    string jName = CrossFluidParam.attrib("species2");

    size_t num = CrossFluidParam.nChildren();
    for (size_t iChild = 0; iChild < num; iChild++) {
        XML_Node& xmlChild = CrossFluidParam.child(iChild);
        string nodeName = toLowerCopy(xmlChild.name());

        if (nodeName == "a_coeff") {
            vector_fp vParams;
            getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");
            string iModel = toLowerCopy(xmlChild.attrib("model"));
            if (iModel == "constant" && vParams.size() == 1) {
                setBinaryCoeffs(iName, jName, vParams[0], 0.0);
            } else if (iModel == "linear_a") {
                setBinaryCoeffs(iName, jName, vParams[0], vParams[1]);
            } else {
                throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid",
                    "unknown model ({}) or wrong number of parameters ({})",
                    iModel, vParams.size());
            }
        }
    }
}

void RedlichKwongMFTP::setParametersFromXML(const XML_Node& thermoNode)
{
    MixtureFugacityTP::setParametersFromXML(thermoNode);
    std::string model = thermoNode["model"];
}

doublereal RedlichKwongMFTP::sresid() const
{
    // note this agrees with tpx
    doublereal rho = density();
    doublereal mmw = meanMolecularWeight();
    doublereal molarV = mmw / rho;
    double hh = m_b_current / molarV;
    doublereal zz = z();
    doublereal dadt = da_dt();
    doublereal T = temperature();
    doublereal sqT = sqrt(T);
    doublereal fac = dadt - m_a_current / (2.0 * T);
    double sresid_mol_R = log(zz*(1.0 - hh)) + log(1.0 + hh) * fac / (sqT * GasConstant * m_b_current);
    return GasConstant * sresid_mol_R;
}

doublereal RedlichKwongMFTP::hresid() const
{
    // note this agrees with tpx
    doublereal rho = density();
    doublereal mmw = meanMolecularWeight();
    doublereal molarV = mmw / rho;
    double hh = m_b_current / molarV;
    doublereal zz = z();
    doublereal dadt = da_dt();
    doublereal T = temperature();
    doublereal sqT = sqrt(T);
    doublereal fac = T * dadt - 3.0 *m_a_current / (2.0);
    return GasConstant * T * (zz - 1.0) + fac * log(1.0 + hh) / (sqT * m_b_current);
}

doublereal RedlichKwongMFTP::liquidVolEst(doublereal TKelvin, doublereal& presGuess) const
{
    double v = m_b_current * 1.1;
    double atmp;
    double btmp;
    calculateAB(TKelvin, atmp, btmp);
    doublereal pres = std::max(psatEst(TKelvin), presGuess);
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

doublereal RedlichKwongMFTP::densityCalc(doublereal TKelvin, doublereal presPa, int phaseRequested, doublereal rhoguess)
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


    doublereal volguess = mmw / rhoguess;
    NSolns_ = solveCubic(TKelvin, presPa, m_a_current, m_b_current, Vroot_);

    doublereal molarVolLast = Vroot_[0];
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

doublereal RedlichKwongMFTP::densSpinodalLiquid() const
{
    double Vroot[3];
    double T = temperature();
    int nsol = solveCubic(T, pressure(), m_a_current, m_b_current, Vroot);
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

    doublereal mmw = meanMolecularWeight();
    return mmw / (0.5 * (vv.first + vv.second));
}

doublereal RedlichKwongMFTP::densSpinodalGas() const
{
    double Vroot[3];
    double T = temperature();
    int nsol = solveCubic(T, pressure(), m_a_current, m_b_current, Vroot);
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

    doublereal mmw = meanMolecularWeight();
    return mmw / (0.5 * (vv.first + vv.second));
}

doublereal RedlichKwongMFTP::dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const
{
    doublereal sqt = sqrt(TKelvin);
    presCalc = GasConstant * TKelvin / (molarVol - m_b_current)
               - m_a_current / (sqt * molarVol * (molarVol + m_b_current));

    doublereal vpb = molarVol + m_b_current;
    doublereal vmb = molarVol - m_b_current;
    doublereal dpdv = (- GasConstant * TKelvin / (vmb * vmb)
                       + m_a_current * (2 * molarVol + m_b_current) / (sqt * molarVol * molarVol * vpb * vpb));
    return dpdv;
}

void RedlichKwongMFTP::pressureDerivatives() const
{
    doublereal TKelvin = temperature();
    doublereal mv = molarVolume();
    doublereal pres;

    dpdV_ = dpdVCalc(TKelvin, mv, pres);
    doublereal sqt = sqrt(TKelvin);
    doublereal vpb = mv + m_b_current;
    doublereal vmb = mv - m_b_current;
    doublereal dadt = da_dt();
    doublereal fac = dadt - m_a_current/(2.0 * TKelvin);
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

void RedlichKwongMFTP::calculateAB(doublereal temp, doublereal& aCalc, doublereal& bCalc) const
{
    bCalc = 0.0;
    aCalc = 0.0;
    if (m_formTempParam == 1) {
        for (size_t i = 0; i < m_kk; i++) {
            bCalc += moleFractions_[i] * b_vec_Curr_[i];
            for (size_t j = 0; j < m_kk; j++) {
                size_t counter = i * m_kk + j;
                doublereal a_vec_Curr = a_coeff_vec(0,counter) + a_coeff_vec(1,counter) * temp;
                aCalc += a_vec_Curr * moleFractions_[i] * moleFractions_[j];
            }
        }
    } else {
        for (size_t i = 0; i < m_kk; i++) {
            bCalc += moleFractions_[i] * b_vec_Curr_[i];
            for (size_t j = 0; j < m_kk; j++) {
                size_t counter = i * m_kk + j;
                doublereal a_vec_Curr = a_coeff_vec(0,counter);
                aCalc += a_vec_Curr * moleFractions_[i] * moleFractions_[j];
            }
        }
    }
}

doublereal RedlichKwongMFTP::da_dt() const
{
    doublereal dadT = 0.0;
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

void RedlichKwongMFTP::calcCriticalConditions(doublereal& pc, doublereal& tc, doublereal& vc) const
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
    doublereal sqrttc, f, dfdt, deltatc;

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

int RedlichKwongMFTP::solveCubic(double T, double pres, double a, double b, double Vroot[3]) const
{

    // Derive the coefficients of the cubic polynomial to solve.
    doublereal an = 1.0;
    doublereal bn = - GasConstant * T / pres;
    doublereal sqt = sqrt(T);
    doublereal cn = - (GasConstant * T * b / pres - a/(pres * sqt) + b * b);
    doublereal dn = - (a * b / (pres * sqt));

    double tmp = a * omega_b / (b * omega_a * GasConstant);
    double pp = 2./3.;
    double tc = pow(tmp, pp);
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    int nSolnValues = MixtureFugacityTP::solveCubic(T, pres, a, b, a, Vroot, an, bn, cn, dn, tc, vc);

    return nSolnValues;
}

}
