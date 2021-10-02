//! @file PengRobinson.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

const double PengRobinson::omega_a = 4.5723552892138218E-01;
const double PengRobinson::omega_b = 7.77960739038885E-02;
const double PengRobinson::omega_vc = 3.07401308698703833E-01;

PengRobinson::PengRobinson(const std::string& infile, const std::string& id_) :
    m_b(0.0),
    m_a(0.0),
    m_aAlpha_mix(0.0),
    m_NSolns(0),
    m_dpdV(0.0),
    m_dpdT(0.0)
{
    fill_n(m_Vroot, 3, 0.0);
    initThermoFile(infile, id_);
}

void PengRobinson::setSpeciesCoeffs(const std::string& species, double a, double b, double w)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
        throw CanteraError("PengRobinson::setSpeciesCoeffs",
            "Unknown species '{}'.", species);
    }

    // Calculate value of kappa (independent of temperature)
    // w is an acentric factor of species
    if (w <= 0.491) {
        m_kappa[k] = 0.37464 + 1.54226*w - 0.26992*w*w;
    } else {
        m_kappa[k] = 0.374642 + 1.487503*w - 0.164423*w*w + 0.016666*w*w*w;
    }

    //Calculate alpha (temperature dependent interaction parameter)
    double critTemp = speciesCritTemperature(a, b); // critical temperature of individual species
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

void PengRobinson::setBinaryCoeffs(const std::string& species_i,
        const std::string& species_j, double a0)
{
    size_t ki = speciesIndex(species_i);
    if (ki == npos) {
        throw CanteraError("PengRobinson::setBinaryCoeffs",
            "Unknown species '{}'.", species_i);
    }
    size_t kj = speciesIndex(species_j);
    if (kj == npos) {
        throw CanteraError("PengRobinson::setBinaryCoeffs",
            "Unknown species '{}'.", species_j);
    }

    m_a_coeffs(ki, kj) = m_a_coeffs(kj, ki) = a0;
    // Calculate alpha_ij
    double alpha_ij = m_alpha[ki] * m_alpha[kj];
    m_aAlpha_binary(ki, kj) = m_aAlpha_binary(kj, ki) = a0*alpha_ij;
}

// ------------Molar Thermodynamic Properties -------------------------

double PengRobinson::cp_mole() const
{
    _updateReferenceStateThermo();
    double T = temperature();
    double mv = molarVolume();
    double vpb = mv + (1 + Sqrt2) * m_b;
    double vmb = mv + (1 - Sqrt2) * m_b;
    calculatePressureDerivatives();
    double cpref = GasConstant * mean_X(m_cp0_R);
    double dHdT_V = cpref + mv * m_dpdT - GasConstant
                    + 1.0 / (2.0 * Sqrt2 * m_b) * log(vpb / vmb) * T * d2aAlpha_dT2();
    return dHdT_V - (mv + T * m_dpdT / m_dpdV) * m_dpdT;
}

double PengRobinson::cv_mole() const
{
    _updateReferenceStateThermo();
    double T = temperature();
    calculatePressureDerivatives();
    return (cp_mole() + T * m_dpdT * m_dpdT / m_dpdV);
}

double PengRobinson::pressure() const
{
    _updateReferenceStateThermo();
    //  Get a copy of the private variables stored in the State object
    double T = temperature();
    double mv = molarVolume();
    double denom = mv * mv + 2 * mv * m_b - m_b * m_b;
    double pp = GasConstant * T / (mv - m_b) - m_aAlpha_mix / denom;
    return pp;
}

double PengRobinson::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void PengRobinson::getActivityCoefficients(double* ac) const
{
    double mv = molarVolume();
    double vpb2 = mv + (1 + Sqrt2) * m_b;
    double vmb2 = mv + (1 - Sqrt2) * m_b;
    double vmb = mv - m_b;
    double pres = pressure();

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
        }
    }
    double num = 0;
    double denom = 2 * Sqrt2 * m_b * m_b;
    double denom2 = m_b * (mv * mv + 2 * mv * m_b - m_b * m_b);
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        num = 2 * m_b * m_pp[k] - m_aAlpha_mix * m_b_coeffs[k];
        ac[k] = (-RTkelvin * log(pres * mv/ RTkelvin) + RTkelvin * log(mv / vmb)
                 + RTkelvin * m_b_coeffs[k] / vmb
                 - (num /denom) * log(vpb2/vmb2)
                 - m_aAlpha_mix * m_b_coeffs[k] * mv/denom2
                );
    }
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(ac[k]/ RTkelvin);
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void PengRobinson::getChemPotentials(double* mu) const
{
    getGibbs_ref(mu);
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RTkelvin * (log(xx));
    }

    double mv = molarVolume();
    double vmb = mv - m_b;
    double vpb2 = mv + (1 + Sqrt2) * m_b;
    double vmb2 = mv + (1 - Sqrt2) * m_b;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
        }
    }
    double pres = pressure();
    double refP = refPressure();
    double denom = 2 * Sqrt2 * m_b * m_b;
    double denom2 = m_b * (mv * mv + 2 * mv * m_b - m_b * m_b);

    for (size_t k = 0; k < m_kk; k++) {
        double num = 2 * m_b * m_pp[k] - m_aAlpha_mix * m_b_coeffs[k];

        mu[k] += (RTkelvin * log(pres/refP) - RTkelvin * log(pres * mv / RTkelvin)
                  + RTkelvin * log(mv / vmb)
                  + RTkelvin * m_b_coeffs[k] / vmb
                  - (num /denom) * log(vpb2/vmb2)
                  - m_aAlpha_mix * m_b_coeffs[k] * mv/denom2
                 );
    }
}

void PengRobinson::getPartialMolarEnthalpies(double* hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar, hbar+m_kk, hbar, RT());
    vector_fp tmp;
    tmp.resize(m_kk,0.0);

    // We calculate m_dpdni
    double T = temperature();
    double mv = molarVolume();
    double vmb = mv - m_b;
    double vpb2 = mv + (1 + Sqrt2) * m_b;
    double vmb2 = mv + (1 - Sqrt2) * m_b;
    double daAlphadT = daAlpha_dT();

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        tmp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            double grad_aAlpha = m_dalphadT[i]/m_alpha[i] + m_dalphadT[k]/m_alpha[k];
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
            tmp[k] +=moleFractions_[i] * m_aAlpha_binary(k, i) * grad_aAlpha;
        }
    }

    double denom = mv * mv + 2 * mv * m_b - m_b * m_b;
    double denom2 = denom * denom;
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        m_dpdni[k] = RTkelvin / vmb + RTkelvin * m_b_coeffs[k] / (vmb * vmb) - 2.0 * m_pp[k] / denom
                    + 2 * vmb * m_aAlpha_mix * m_b_coeffs[k] / denom2;
    }

    double fac = T * daAlphadT - m_aAlpha_mix;
    calculatePressureDerivatives();
    double fac2 = mv + T * m_dpdT / m_dpdV;
    double fac3 = 2 * Sqrt2 * m_b * m_b;
    double fac4 = 0;
    for (size_t k = 0; k < m_kk; k++) {
        fac4 = T*tmp[k] -2 * m_pp[k];
        double hE_v = mv * m_dpdni[k] - RTkelvin - m_b_coeffs[k] / fac3  * log(vpb2 / vmb2) * fac
                     + (mv * m_b_coeffs[k]) /(m_b * denom) * fac
                     + 1/(2 * Sqrt2 * m_b) * log(vpb2 / vmb2) * fac4;
        hbar[k] = hbar[k] + hE_v;
        hbar[k] -= fac2 * m_dpdni[k];
    }
}

void PengRobinson::getPartialMolarEntropies(double* sbar) const
{
    // Using the identity : (hk - T*sk) = gk
    double T = temperature();
    getPartialMolarEnthalpies(sbar);
    getChemPotentials(m_tmpV.data());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] = (sbar[k] - m_tmpV[k])/T;
    }
}

void PengRobinson::getPartialMolarIntEnergies(double* ubar) const
{
    // u_i = h_i - p*v_i
    double p = pressure();
    getPartialMolarEnthalpies(ubar);
    getPartialMolarVolumes(m_tmpV.data());
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = ubar[k] - p*m_tmpV[k];
    }
}

void PengRobinson::getPartialMolarCp(double* cpbar) const
{
    throw NotImplementedError("PengRobinson::getPartialMolarCp");
}

void PengRobinson::getPartialMolarVolumes(double* vbar) const
{
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
        }
    }

    double mv = molarVolume();
    double vmb = mv - m_b;
    double vpb = mv + m_b;
    double fac = mv * mv + 2 * mv * m_b - m_b * m_b;
    double fac2 = fac * fac;
    double RTkelvin = RT();

    for (size_t k = 0; k < m_kk; k++) {
        double num = (RTkelvin + RTkelvin * m_b/ vmb + RTkelvin * m_b_coeffs[k] / vmb
                      + RTkelvin * m_b * m_b_coeffs[k] /(vmb * vmb)
                      - 2 * mv * m_pp[k] / fac
                      + 2 * mv * vmb * m_aAlpha_mix * m_b_coeffs[k] / fac2
                     );
        double denom = (pressure() + RTkelvin * m_b / (vmb * vmb)
                        + m_aAlpha_mix/fac
                        - 2 * mv* vpb * m_aAlpha_mix / fac2
                       );
        vbar[k] = num / denom;
    }
}

double PengRobinson::speciesCritTemperature(double a, double b) const
{
    if (b <= 0.0) {
        return 1000000.;
    } else if (a <= 0.0) {
        return 0.0;
    } else {
        return a * omega_b / (b * omega_a * GasConstant);
    }
}

bool PengRobinson::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        m_a_coeffs.resize(m_kk, m_kk, 0.0);
        m_b_coeffs.push_back(0.0);
        m_aAlpha_binary.resize(m_kk, m_kk, 0.0);
        m_kappa.push_back(0.0);

        m_alpha.push_back(0.0);
        m_dalphadT.push_back(0.0);
        m_d2alphadT2.push_back(0.0);

        m_pp.push_back(0.0);
        m_partialMolarVolumes.push_back(0.0);
        m_dpdni.push_back(0.0);
    }
    return added;
}

vector<double> PengRobinson::getCoeff(const std::string& iName)
{
    vector_fp spCoeff{ NAN, NAN, NAN };

    // Get number of species in the database
    // open xml file critProperties.xml
    XML_Node* doc = get_XML_File("critProperties.xml");
    size_t nDatabase = doc->nChildren();

    // Loop through all species in the database and attempt to match supplied
    // species to each. If present, calculate pureFluidParameters a_k and b_k
    // based on crit properties T_c and P_c:
    for (size_t isp = 0; isp < nDatabase; isp++) {
        XML_Node& acNodeDoc = doc->child(isp);
        std::string iNameLower = toLowerCopy(iName);
        std::string dbName = toLowerCopy(acNodeDoc.attrib("name"));

        // Attempt to match provided species iName to current database species
        //  dbName:
        if (iNameLower == dbName) {
            // Read from database and calculate a and b coefficients
            double vParams;
            double T_crit = 0.0, P_crit = 0.0, w_ac = 0.0;

            if (acNodeDoc.hasChild("Tc")) {
                vParams = 0.0;
                XML_Node& xmlChildCoeff = acNodeDoc.child("Tc");
                if (xmlChildCoeff.hasAttrib("value")) {
                    std::string critTemp = xmlChildCoeff.attrib("value");
                    vParams = strSItoDbl(critTemp);
                }
                if (vParams <= 0.0) { //Assuming that Pc and Tc are non zero.
                    throw CanteraError("PengRobinson::getCoeff",
                        "Critical Temperature must be positive");
                }
                T_crit = vParams;
            }
            if (acNodeDoc.hasChild("Pc")) {
                vParams = 0.0;
                XML_Node& xmlChildCoeff = acNodeDoc.child("Pc");
                if (xmlChildCoeff.hasAttrib("value")) {
                    std::string critPressure = xmlChildCoeff.attrib("value");
                    vParams = strSItoDbl(critPressure);
                }
                if (vParams <= 0.0) { //Assuming that Pc and Tc are non zero.
                    throw CanteraError("PengRobinson::getCoeff",
                        "Critical Pressure must be positive");
                }
                P_crit = vParams;
            }
            if (acNodeDoc.hasChild("omega")) {
                vParams = 0.0;
                XML_Node& xmlChildCoeff = acNodeDoc.child("omega");
                if (xmlChildCoeff.hasChild("value")) {
                    std::string acentric_factor = xmlChildCoeff.attrib("value");
                    vParams = strSItoDbl(acentric_factor);
                }
                w_ac = vParams;
            }

            spCoeff[0] = omega_a * (GasConstant * GasConstant) * (T_crit * T_crit) / P_crit; //coeff a
            spCoeff[1] = omega_b * GasConstant * T_crit / P_crit; // coeff b
            spCoeff[2] = w_ac; // acentric factor
            break;
        }
    }
    // If the species is not present in the database, throw an error
    if(isnan(spCoeff[0]))
    {
        throw CanteraError("PengRobinson::getCoeff",
            "Species '{}' is not present in the database", iName);
    }
    return spCoeff;
}

void PengRobinson::initThermo()
{
    for (auto& item : m_species) {
        // Read a and b coefficients and acentric factor w_ac from species input
        // information, specified in a YAML input file.
        if (item.second->input.hasKey("equation-of-state")) {
            auto eos = item.second->input["equation-of-state"].getMapWhere(
                "model", "Peng-Robinson");
            double a0 = 0;
            if (eos["a"].isScalar()) {
                a0 = eos.convert("a", "Pa*m^6/kmol^2");
            }
            double b = eos.convert("b", "m^3/kmol");
            // unitless acentric factor:
            double w = eos["acentric-factor"].asDouble();

            setSpeciesCoeffs(item.first, a0, b, w);
            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& item2 : binary_a) {
                    double a0 = 0;
                    if (item2.second.isScalar()) {
                        a0 = units.convert(item2.second, "Pa*m^6/kmol^2");
                    }
                    setBinaryCoeffs(item.first, item2.first, a0);
                }
            }
        } else {
            // Check if a and b are already populated for this species (only the
            // diagonal elements of a). If not, then search 'critProperties.xml'
            // to find critical temperature and pressure to calculate a and b.
            size_t k = speciesIndex(item.first);
            if (m_a_coeffs(k, k) == 0.0) {
                vector<double> coeffs = getCoeff(item.first);

                // Check if species was found in the database of critical
                // properties, and assign the results
                if (!isnan(coeffs[0])) {
                    setSpeciesCoeffs(item.first, coeffs[0], coeffs[1], coeffs[2]);
                }
            }
        }
    }
}

double PengRobinson::sresid() const
{
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

double PengRobinson::hresid() const
{
    double molarV = molarVolume();
    double zz = z();
    double aAlpha_1 = daAlpha_dT();
    double T = temperature();
    double vpb = molarV + (1 + Sqrt2) * m_b;
    double vmb = molarV + (1 - Sqrt2) * m_b;
    double fac = 1 / (2.0 * Sqrt2 * m_b);
    return GasConstant * T * (zz - 1.0) + fac * log(vpb / vmb) * (T * aAlpha_1 - m_aAlpha_mix);
}

double PengRobinson::liquidVolEst(double T, double& presGuess) const
{
    double v = m_b * 1.1;
    double atmp;
    double btmp;
    double aAlphatmp;
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

double PengRobinson::densityCalc(double T, double presPa, int phaseRequested, double rhoGuess)
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
        if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT || phaseRequested == FLUID_UNDEFINED) {
            molarVolLast = m_Vroot[0];
        } else {
            return -2.0;
        }
    } else if (m_NSolns == -1) {
        if (phaseRequested >= FLUID_LIQUID_0 || phaseRequested == FLUID_UNDEFINED || phaseRequested == FLUID_SUPERCRIT) {
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

double PengRobinson::densSpinodalLiquid() const
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

double PengRobinson::densSpinodalGas() const
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

double PengRobinson::dpdVCalc(double T, double molarVol, double& presCalc) const
{
    double denom = molarVol * molarVol + 2 * molarVol * m_b - m_b * m_b;
    double vpb = molarVol + m_b;
    double vmb = molarVol - m_b;
    double dpdv = -GasConstant * T / (vmb * vmb) + 2 * m_aAlpha_mix * vpb / (denom*denom);
    return dpdv;
}

void PengRobinson::calculatePressureDerivatives() const
{
    double T = temperature();
    double mv = molarVolume();
    double pres;

    m_dpdV = dpdVCalc(T, mv, pres);
    double vmb = mv - m_b;
    double denom = mv * mv + 2 * mv * m_b - m_b * m_b;
    m_dpdT = (GasConstant / vmb - daAlpha_dT() / denom);
}

void PengRobinson::updateMixingExpressions()
{
    double temp = temperature();

    // Update individual alpha
    for (size_t j = 0; j < m_kk; j++) {
        double critTemp_j = speciesCritTemperature(m_a_coeffs(j,j), m_b_coeffs[j]);
        double sqt_alpha = 1 + m_kappa[j] * (1 - sqrt(temp / critTemp_j));
        m_alpha[j] = sqt_alpha*sqt_alpha;
    }

    //Update aAlpha_i, j
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            m_aAlpha_binary(i, j) = sqrt(m_alpha[i] * m_alpha[j]) * m_a_coeffs(i,j);
        }
    }
    calculateAB(m_a,m_b,m_aAlpha_mix);
}

void PengRobinson::calculateAB(double& aCalc, double& bCalc, double& aAlphaCalc) const
{
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

double PengRobinson::daAlpha_dT() const
{
    double daAlphadT = 0.0, k, Tc, sqtTr, coeff1, coeff2;
    for (size_t i = 0; i < m_kk; i++) {
        // Calculate first derivative of alpha for individual species
        Tc = speciesCritTemperature(m_a_coeffs(i,i), m_b_coeffs[i]);
        sqtTr = sqrt(temperature() / Tc); //we need species critical temperature
        coeff1 = 1 / (Tc*sqtTr);
        coeff2 = sqtTr - 1;
        k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
    }
    //Calculate mixture derivative
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            daAlphadT += moleFractions_[i] * moleFractions_[j] * 0.5 * m_aAlpha_binary(i, j)
                                             * (m_dalphadT[i] / m_alpha[i] + m_dalphadT[j] / m_alpha[j]);
        }
    }
    return daAlphadT;
}

double PengRobinson::d2aAlpha_dT2() const
{
    for (size_t i = 0; i < m_kk; i++) {
        double Tcrit_i = speciesCritTemperature(m_a_coeffs(i, i), m_b_coeffs[i]);
        double sqt_Tr = sqrt(temperature() / Tcrit_i); //we need species critical temperature
        double coeff1 = 1 / (Tcrit_i*sqt_Tr);
        double coeff2 = sqt_Tr - 1;
        //  Calculate first and second derivatives of alpha for individual species
        double k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
        m_d2alphadT2[i] = (k*k + k) * coeff1 / (2*sqt_Tr*sqt_Tr*Tcrit_i);
    }

    //Calculate mixture derivative
    double d2aAlphadT2 = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        double alphai = m_alpha[i];
        for (size_t j = 0; j < m_kk; j++) {
            double alphaj = m_alpha[j];
            double alphaij = alphai * alphaj;
            double term1 = m_d2alphadT2[i] / alphai + m_d2alphadT2[j] / alphaj;
            double term2 = 2 * m_dalphadT[i] * m_dalphadT[j] / alphaij;
            double term3 = m_dalphadT[i] / alphai + m_dalphadT[j] / alphaj;
            d2aAlphadT2 += 0.5 * moleFractions_[i] * moleFractions_[j] * m_aAlpha_binary(i, j)
                                       * (term1 + term2 - 0.5 * term3 * term3);
        }
    }
    return d2aAlphadT2;
}

void PengRobinson::calcCriticalConditions(double& pc, double& tc, double& vc) const
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

int PengRobinson::solveCubic(double T, double pres, double a, double b, double aAlpha, double Vroot[3]) const
{
    // Derive the coefficients of the cubic polynomial (in terms of molar volume v) to solve.
    double bsqr = b * b;
    double RT_p = GasConstant * T / pres;
    double aAlpha_p = aAlpha / pres;
    double an = 1.0;
    double bn = (b - RT_p);
    double cn = -(2 * RT_p * b - aAlpha_p + 3 * bsqr);
    double dn = (bsqr * RT_p + bsqr * b - aAlpha_p * b);

    double tc = a * omega_b / (b * omega_a * GasConstant);
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    int nSolnValues = MixtureFugacityTP::solveCubic(T, pres, a, b, aAlpha, Vroot, an, bn, cn, dn, tc, vc);

    return nSolnValues;
}

}
