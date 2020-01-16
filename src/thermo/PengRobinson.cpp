//! @file PengRobinson.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <boost/math/tools/roots.hpp>

#define _USE_MATH_DEFINES
#include<math.h>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

const double PengRobinson::omega_a = 4.5723552892138218E-01;
const double PengRobinson::omega_b = 7.77960739038885E-02;
const double PengRobinson::omega_vc = 3.07401308698703833E-01;

PengRobinson::PengRobinson() :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    m_aAlpha_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    fill_n(Vroot_, 3, 0.0);
}

PengRobinson::PengRobinson(const std::string& infile, const std::string& id_) :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    m_aAlpha_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    fill_n(Vroot_, 3, 0.0);
    initThermoFile(infile, id_);
}

/*PengRobinson::PengRobinson(XML_Node& phaseRefRoot, const std::string& id_) :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    m_aAlpha_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    fill_n(Vroot_, 3, 0.0);
    importPhase(phaseRefRoot, this);
}
*/
void PengRobinson::calculateAlpha(const std::string& species, double a, double b, double w)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
        throw CanteraError("PengRobinson::calculateAlpha",
            "Unknown species '{}'.", species);
    }

    // Calculate value of kappa (independent of temperature)
    // w is an acentric factor of species and must be specified in the CTI file
    if (isnan(w)){
        throw CanteraError("PengRobinson::calculateALpha",
            "No acentric factor loaded.");
    } else if (w <= 0.491) {
        kappa_vec_[k] = 0.37464 + 1.54226*w - 0.26992*w*w;
    } else {
        kappa_vec_[k] = 0.374642 + 1.487503*w - 0.164423*w*w + 0.016666*w*w*w;
    }

    //Calculate alpha (temperature dependent interaction parameter)
    double critTemp = speciesCritTemperature(a, b); // critical temperature of individual species
    double sqt_T_r = sqrt(temperature() / critTemp);
    double sqt_alpha = 1 + kappa_vec_[k] * (1 - sqt_T_r);
    alpha_vec_Curr_[k] = sqt_alpha*sqt_alpha;
}

void PengRobinson::setSpeciesCoeffs(const std::string& species,
                                        double a, double b, double w)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
        throw CanteraError("PengRobinson::setSpeciesCoeffs",
            "Unknown species '{}'.", species);
    }
    size_t counter = k + m_kk * k;
    a_coeff_vec(0, counter) = a;
    // we store this locally because it is used below to calculate a_Alpha:
    double aAlpha_k = a*alpha_vec_Curr_[k];
    aAlpha_coeff_vec(0, counter) = aAlpha_k;

    // standard mixing rule for cross-species interaction term
    for (size_t j = 0; j < m_kk; j++) {
        if (k == j) {
            continue;
        }
        double a0kj = sqrt(a_coeff_vec(0, j + m_kk * j) * a);
        double aAlpha_j = a*alpha_vec_Curr_[j];
        double a_Alpha = sqrt(aAlpha_j*aAlpha_k);
        if (a_coeff_vec(0, j + m_kk * k) == 0) {
            a_coeff_vec(0, j + m_kk * k) = a0kj;
            aAlpha_coeff_vec(0, j + m_kk * k) = a_Alpha;
            a_coeff_vec(0, k + m_kk * j) = a0kj;
            aAlpha_coeff_vec(0, k + m_kk * j) = a_Alpha;
        }
    }
    a_coeff_vec.getRow(0, a_vec_Curr_.data());
    aAlpha_coeff_vec.getRow(0, aAlpha_vec_Curr_.data());
    b_vec_Curr_[k] = b;
}

void PengRobinson::setBinaryCoeffs(const std::string& species_i,
        const std::string& species_j, double a0, double alpha)
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

    size_t counter1 = ki + m_kk * kj;
    size_t counter2 = kj + m_kk * ki;
    a_coeff_vec(0, counter1) = a_coeff_vec(0, counter2) = a0;
    aAlpha_coeff_vec(0, counter1) = aAlpha_coeff_vec(0, counter2) = a0*alpha;
    a_vec_Curr_[counter1] = a_vec_Curr_[counter2] = a0;
    aAlpha_vec_Curr_[counter1] = aAlpha_vec_Curr_[counter2] = a0*alpha;
}

// ------------Molar Thermodynamic Properties -------------------------

double PengRobinson::cp_mole() const
{
    _updateReferenceStateThermo();
    double T = temperature();
    double mv = molarVolume();
    double vpb = mv + (1 + M_SQRT2)*m_b_current;
    double vmb = mv + (1 - M_SQRT2)*m_b_current;
    pressureDerivatives();
    double cpref = GasConstant * mean_X(m_cp0_R);
    double dHdT_V = cpref + mv * dpdT_ - GasConstant
                        + 1.0 / (2.0 * M_SQRT2 *m_b_current) * log(vpb / vmb) * T *d2aAlpha_dT2();
    return dHdT_V - (mv + T * dpdT_ / dpdV_) * dpdT_;
}

double PengRobinson::cv_mole() const
{
    _updateReferenceStateThermo();
    double T = temperature();
    pressureDerivatives();
    return (cp_mole() + T* dpdT_* dpdT_ / dpdV_);
}

double PengRobinson::pressure() const
{
    _updateReferenceStateThermo();
    //  Get a copy of the private variables stored in the State object
    double T = temperature();
    double mv = molarVolume();
    double denom = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    double pp = GasConstant * T / (mv - m_b_current) - m_aAlpha_current / denom;
    return pp;
}

void PengRobinson::setTemperature(const double temp)
{
    Phase::setTemperature(temp);
    _updateReferenceStateThermo();
    updateAB();
}

void PengRobinson::compositionChanged()
{
    MixtureFugacityTP::compositionChanged();
    updateAB();
}

double PengRobinson::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void PengRobinson::getActivityCoefficients(double* ac) const
{
    double mv = molarVolume();
    //double T = temperature();
    double vpb2 = mv + (1 + M_SQRT2)*m_b_current;
    double vmb2 = mv + (1 - M_SQRT2)*m_b_current;
    double vmb = mv - m_b_current;
    double pres = pressure();

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }
    double num = 0;
    double den = 2 * M_SQRT2 * m_b_current * m_b_current;
    double den2 = m_b_current*(mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current);
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        num = 2 * m_b_current * m_pp[k] - m_aAlpha_current* b_vec_Curr_[k];
        ac[k] = (-RTkelvin *log(pres*mv/ RTkelvin) + RTkelvin * log(mv / vmb)
                 + RTkelvin * b_vec_Curr_[k] / vmb
                 - (num /den) * log(vpb2/vmb2)
                 - m_aAlpha_current* b_vec_Curr_[k] * mv/den2
                );
    }
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(ac[k]/ RTkelvin);
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void PengRobinson::getChemPotentials_RT(double* muRT) const
{
    getChemPotentials(muRT);
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= 1.0 / RTkelvin;
    }
}

void PengRobinson::getChemPotentials(double* mu) const
{
    getGibbs_ref(mu);
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RTkelvin *(log(xx));
    }

    double mv = molarVolume();
    double vmb = mv - m_b_current;
    double vpb2 = mv + (1 + M_SQRT2)*m_b_current;
    double vmb2 = mv + (1 - M_SQRT2)*m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }
    double pres = pressure();
    double refP = refPressure();
    double num = 0;
    double den = 2 * M_SQRT2 * m_b_current * m_b_current;
    double den2 = m_b_current*(mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current);

    for (size_t k = 0; k < m_kk; k++) {
        num = 2 * m_b_current * m_pp[k] - m_aAlpha_current* b_vec_Curr_[k];

        mu[k] += (RTkelvin * log(pres/refP) - RTkelvin * log(pres * mv / RTkelvin)
                  + RTkelvin * log(mv / vmb)
                  + RTkelvin * b_vec_Curr_[k] / vmb
                  - (num /den) * log(vpb2/vmb2)
                  - m_aAlpha_current* b_vec_Curr_[k] * mv/den2
                 );
    }
}

void PengRobinson::getPartialMolarEnthalpies(double* hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar, hbar+m_kk, hbar, RT());

    // We calculate dpdni_
    double T = temperature();
    double mv = molarVolume();
    double vmb = mv - m_b_current;
    double vpb2 = mv + (1 + M_SQRT2)*m_b_current;
    double vmb2 = mv + (1 - M_SQRT2)*m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }

    double den = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    double den2 = den*den;
    double RTkelvin = RT();
    for (size_t k = 0; k < m_kk; k++) {
        dpdni_[k] = RTkelvin /vmb + RTkelvin * b_vec_Curr_[k] / (vmb * vmb) - 2.0 * m_pp[k] / den
                    + 2 * vmb * m_aAlpha_current * b_vec_Curr_[k] / den2;
    }

    double daAlphadT = daAlpha_dT();
    double fac = T * daAlphadT - m_aAlpha_current;

    pressureDerivatives();
    double fac2 = mv + T * dpdT_ / dpdV_;
        double fac3 = 2 * M_SQRT2 * m_b_current *m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
        double hE_v = mv * dpdni_[k] - RTkelvin + (2 * m_b_current - b_vec_Curr_[k]) / fac3  * log(vpb2 / vmb2)*fac
                    + (mv * b_vec_Curr_[k]) /(m_b_current*den) * fac;
        hbar[k] = hbar[k] + hE_v;
        hbar[k] -= fac2 * dpdni_[k];
    }
}

void PengRobinson::getPartialMolarEntropies(double* sbar) const
{
    getEntropy_R_ref(sbar);
    scale(sbar, sbar+m_kk, sbar, GasConstant);
    double T = temperature();
    double mv = molarVolume();
    double vmb = mv - m_b_current;
    double vpb2 = mv + (1 + M_SQRT2)*m_b_current;
    double vmb2 = mv + (1 - M_SQRT2)*m_b_current;
    double refP = refPressure();
    double daAlphadT = daAlpha_dT();
    double coeff1 = 0;
    double den1 = 2 * M_SQRT2 * m_b_current * m_b_current;
    double den2 = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;

    // Calculate sum(n_j (a alpha)_i, k * (1/alpha_k d/dT(alpha_k))) -> m_pp
    // Calculate sum(n_j (a alpha)_i, k * (1/alpha_i d/dT(alpha_i))) -> m_tmpV
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        m_tmpV[k] = 0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
                        m_tmpV[k] += moleFractions_[i] * a_coeff_vec(1, counter) *(dalphadT_vec_Curr_[i] / alpha_vec_Curr_[i]);
        }
        m_pp[k] = m_pp[k] * dalphadT_vec_Curr_[k] / alpha_vec_Curr_[k];
    }


    for (size_t k = 0; k < m_kk; k++) {
        coeff1 = m_b_current * (m_pp[k] + m_tmpV[k]) - daAlphadT * b_vec_Curr_[k];
        sbar[k] += GasConstant * log(GasConstant * T / (refP * mv))
                + GasConstant
                + GasConstant * log(mv / vmb)
                + GasConstant * b_vec_Curr_[k] / vmb
                - coeff1* log(vpb2 / vmb2) / den1
                - b_vec_Curr_[k] * mv * daAlphadT / den2 / m_b_current;
    }
    pressureDerivatives();
    getPartialMolarVolumes(m_partialMolarVolumes.data());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= m_partialMolarVolumes[k] * dpdT_;
    }
}

void PengRobinson::getPartialMolarIntEnergies(double* ubar) const
{
    getIntEnergy_RT(ubar);
    scale(ubar, ubar+m_kk, ubar, RT());
}

void PengRobinson::getPartialMolarCp(double* cpbar) const
{
    getCp_R(cpbar);
    scale(cpbar, cpbar+m_kk, cpbar, GasConstant);
}

void PengRobinson::getPartialMolarVolumes(double* vbar) const
{
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }

    double mv = molarVolume();
    double vmb = mv - m_b_current;
    double vpb = mv + m_b_current;
    double fac = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    double fac2 = fac * fac;
    double RTkelvin = RT(); 

    for (size_t k = 0; k < m_kk; k++) {
        double num = (RTkelvin + RTkelvin * m_b_current/ vmb + RTkelvin * b_vec_Curr_[k] / vmb
                          + RTkelvin * m_b_current * b_vec_Curr_[k] /(vmb * vmb)
                          - 2 * mv * m_pp[k] / fac
                          + 2 * mv * vmb * m_aAlpha_current * b_vec_Curr_[k] / fac2
                         );
        double denom = (pressure() + RTkelvin * m_b_current / (vmb * vmb)
                            + m_aAlpha_current/fac
                            - 2 * mv* vpb *m_aAlpha_current / fac2
                         );
        vbar[k] = num / denom;
    }
}

double PengRobinson::speciesCritTemperature(double a, double b) const
{
    double pc, tc, vc;
    calcCriticalConditions(a, b, pc, tc, vc);
    return tc;
}

double PengRobinson::critTemperature() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return tc;
}

double PengRobinson::critPressure() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return pc;
}

double PengRobinson::critVolume() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return vc;
}

double PengRobinson::critCompressibility() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return pc*vc/tc/GasConstant;
}

double PengRobinson::critDensity() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    double mmw = meanMolecularWeight();
    return mmw / vc;
}

void PengRobinson::setToEquilState(const double* mu_RT)
{
    double tmp, tmp2;
    _updateReferenceStateThermo();
    getGibbs_RT_ref(m_tmpV.data());

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    double pres = 0.0;
    double m_p0 = refPressure();
    for (size_t k = 0; k < m_kk; k++) {
        tmp = -m_tmpV[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setState_PX(pres, &m_pp[0]);
}

bool PengRobinson::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        a_vec_Curr_.resize(m_kk * m_kk, 0.0);
        b_vec_Curr_.push_back(0.0);
        a_vec_Curr_.push_back(0.0);
        aAlpha_vec_Curr_.resize(m_kk * m_kk, 0.0);
        aAlpha_vec_Curr_.push_back(0.0);
        kappa_vec_.push_back(0.0);

        alpha_vec_Curr_.push_back(0.0);
        a_coeff_vec.resize(1, m_kk * m_kk, 0.0);
        aAlpha_coeff_vec.resize(1, m_kk * m_kk, 0.0);
        dalphadT_vec_Curr_.push_back(0.0);
        d2alphadT2_.push_back(0.0);

        m_pp.push_back(0.0);
        m_tmpV.push_back(0.0);
        m_partialMolarVolumes.push_back(0.0);
        dpdni_.push_back(0.0);
    }
    return added;
}

vector<double> PengRobinson::getCoeff(const std::string& iName)
{
        vector_fp spCoeff{ NAN, NAN };

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
                double T_crit, P_crit = 0.0;

                if (acNodeDoc.hasChild("Tc")) {
                    vParams = 0.0;
                    XML_Node& xmlChildCoeff = acNodeDoc.child("Tc");
                    if (xmlChildCoeff.hasAttrib("value")) {
                        std::string critTemp = xmlChildCoeff.attrib("value");
                        vParams = strSItoDbl(critTemp);
                    }
                    if (vParams <= 0.0) { //Assuming that Pc and Tc are non zero.
                        throw CanteraError("PengRobinson::getCoeff",
                            "Critical Temperature must be positive ");
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
                            "Critical Pressure must be positive ");
                    }
                    P_crit = vParams;
                }

                //Assuming no temperature dependence
                spCoeff[0] = omega_a * (GasConstant* GasConstant) * (T_crit* T_crit) / P_crit; //coeff a
                spCoeff[1] = omega_b * GasConstant * T_crit / P_crit; // coeff b
                break;
            }
        }
        return spCoeff;
}

void PengRobinson::initThermo()
{
    for (auto& item : m_species) {
        // Read a and b coefficients and acentric factor w_ac from species 'input'
        // information (i.e. as specified in a YAML input file)
        if (item.second->input.hasKey("equation-of-state")) {
            auto eos = item.second->input["equation-of-state"].as<AnyMap>();
            if (eos.getString("model", "") != "Peng-Robinson") {
                throw InputFileError("PengRobinson::initThermo", eos,
                    "Expected species equation of state to be 'Peng-Robinson', "
                    "but got '{}' instead", eos.getString("model", ""));
            }
            double a0 = 0, a1 = 0;
            if (eos["a"].isScalar()) {
                a0 = eos.convert("a", "Pa*m^6/kmol^2");
            } else {
                auto avec = eos["a"].asVector<AnyValue>(2);
                a0 = eos.units().convert(avec[0], "Pa*m^6/kmol^2");
                a1 = eos.units().convert(avec[1], "Pa*m^6/kmol^2/K");
            }
            double b = eos.convert("b", "m^3/kmol");
            // unitless acentric factor:
            double w = eos.getDouble("w_ac",NAN);

            calculateAlpha(item.first, a0, b, w);
            setSpeciesCoeffs(item.first, a0, b, w);
            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& item2 : binary_a) {
                    double a0 = 0, a1 = 0;
                    if (item2.second.isScalar()) {
                        a0 = units.convert(item2.second, "Pa*m^6/kmol^2");
                    } else {
                        auto avec = item2.second.asVector<AnyValue>(2);
                        a0 = units.convert(avec[0], "Pa*m^6/kmol^2");
                        a1 = units.convert(avec[1], "Pa*m^6/kmol^2/K");
                    }
                    setBinaryCoeffs(item.first, item2.first, a0, a1);
                }
            }
        } else {
            // Check if a and b are already populated for this species (only the
            // diagonal elements of a). If not, then search 'critProperties.xml'
            // to find critical temperature and pressure to calculate a and b.
            size_t k = speciesIndex(item.first);
            if (isnan(a_coeff_vec(0, k + m_kk * k))) {
                // coeffs[0] = a0, coeffs[1] = b;
                vector<double> coeffs = getCoeff(item.first);

                // Check if species was found in the database of critical
                // properties, and assign the results
                if (!isnan(coeffs[0])) {
                    // Assuming no temperature dependence (i.e. a1 = 0)
                    setSpeciesCoeffs(item.first, coeffs[0], 0.0, coeffs[1]);
                }
            }
        }
    }
}

double PengRobinson::sresid() const
{
    double molarV = molarVolume();
    double hh = m_b_current / molarV;
    double zz = z();
    double alpha_1 = daAlpha_dT();
    double vpb = molarV + (1.0 + M_SQRT2) *m_b_current;
    double vmb = molarV + (1.0 - M_SQRT2) *m_b_current;
    double fac = alpha_1 / (2.0 * M_SQRT2 * m_b_current);
    double sresid_mol_R = log(zz*(1.0 - hh)) + fac * log(vpb / vmb) / GasConstant;
    return GasConstant * sresid_mol_R;
}

double PengRobinson::hresid() const
{
    double molarV = molarVolume();
    double zz = z();
    double aAlpha_1 = daAlpha_dT();
    double T = temperature();
    double vpb = molarV + (1 + M_SQRT2) *m_b_current;
    double vmb = molarV + (1 - M_SQRT2) *m_b_current;
    double fac = 1 / (2.0 * M_SQRT2 * m_b_current);
    return GasConstant * T * (zz - 1.0) + fac * log(vpb / vmb) *(T * aAlpha_1 - m_aAlpha_current);
}

double PengRobinson::liquidVolEst(double T, double& presGuess) const
{
    double v = m_b_current * 1.1;
    double atmp;
    double btmp;
    double aAlphatmp;
    calculateAB(T, atmp, btmp, aAlphatmp);
    double pres = std::max(psatEst(T), presGuess);
    double Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
                int nsol = NicholsCall(T, pres, atmp, btmp, aAlphatmp, Vroot);
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
    // It's necessary to set the temperature so that m_aAlpha_current is set correctly.
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
    NSolns_ = NicholsCall(T, presPa, m_a_current, m_b_current, m_aAlpha_current, Vroot_);

    double molarVolLast = Vroot_[0];
    if (NSolns_ >= 2) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            molarVolLast = Vroot_[0];
        } else if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[2];
        } else {
            if (volGuess > Vroot_[1]) {
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
        } else if (T > tcrit) {
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

double PengRobinson::densSpinodalLiquid() const
{
    double Vroot[3];
    double T = temperature();
    int nsol = NicholsCall(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
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
    int nsol = NicholsCall(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
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

double PengRobinson::pressureCalc(double T, double molarVol) const
{
    double den = molarVol * molarVol + 2 * molarVol * m_b_current - m_b_current * m_b_current;
    double pres = GasConstant * T / (molarVol - m_b_current) - m_aAlpha_current / den;
    return pres;
}

double PengRobinson::dpdVCalc(double T, double molarVol, double& presCalc) const
{
    double den = molarVol * molarVol + 2 * molarVol * m_b_current - m_b_current * m_b_current;
    presCalc = GasConstant * T / (molarVol - m_b_current) - m_aAlpha_current/ den;

    double vpb = molarVol + m_b_current;
    double vmb = molarVol - m_b_current;
    double dpdv = -GasConstant * T / (vmb * vmb) + 2 *m_aAlpha_current * vpb / (den*den);
    return dpdv;
}

void PengRobinson::pressureDerivatives() const
{
    double T = temperature();
    double mv = molarVolume();
    double pres;

    dpdV_ = dpdVCalc(T, mv, pres);
    double vmb = mv - m_b_current;
    double den = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    dpdT_ = (GasConstant / vmb - daAlpha_dT() / den);
}

void PengRobinson::updateMixingExpressions()
{
    updateAB();
}

void PengRobinson::updateAB()
{
    double temp = temperature();
    //Update aAlpha_i
    double sqt_alpha;

    // Update indiviual alpha
    for (size_t j = 0; j < m_kk; j++) {
        size_t counter = j * m_kk + j;
        double critTemp_j = speciesCritTemperature(a_vec_Curr_[counter],b_vec_Curr_[j]);
        sqt_alpha = 1 + kappa_vec_[j] * (1 - sqrt(temp / critTemp_j));
        alpha_vec_Curr_[j] = sqt_alpha*sqt_alpha;
    }

    //Update aAlpha_i, j
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            size_t counter = i * m_kk + j;
            a_vec_Curr_[counter] = a_coeff_vec(0, counter);
            aAlpha_vec_Curr_[counter] = sqrt(alpha_vec_Curr_[i] * alpha_vec_Curr_[j]) * a_coeff_vec(0, counter);
        }
    }

    m_b_current = 0.0;
    m_a_current = 0.0;
    m_aAlpha_current = 0.0;

    for (size_t i = 0; i < m_kk; i++) {
        m_b_current += moleFractions_[i] * b_vec_Curr_[i];
        for (size_t j = 0; j < m_kk; j++) {
            m_a_current += a_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
                        m_aAlpha_current += aAlpha_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
        }
    }
}

void PengRobinson::calculateAB(double temp, double& aCalc, double& bCalc, double& aAlphaCalc) const
{
    bCalc = 0.0;
    aCalc = 0.0;
    aAlphaCalc = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        bCalc += moleFractions_[i] * b_vec_Curr_[i];
        for (size_t j = 0; j < m_kk; j++) {
            size_t counter = i * m_kk + j;
            double a_vec_Curr = a_coeff_vec(0, counter);
            aCalc += a_vec_Curr * moleFractions_[i] * moleFractions_[j];
            aAlphaCalc += aAlpha_vec_Curr_[counter] * moleFractions_[i] * moleFractions_[j];
        }
    }
}

double PengRobinson::daAlpha_dT() const
{
    double daAlphadT = 0.0, temp, k, Tc = 0.0, sqtTr = 0.0;
    double coeff1, coeff2;
    for (size_t i = 0; i < m_kk; i++) {
        size_t counter = i + m_kk * i;
        // Calculate first derivative of alpha for individual species
        Tc = speciesCritTemperature(a_vec_Curr_[counter], b_vec_Curr_[i]);
        sqtTr = sqrt(temperature() / Tc); //we need species critical temperature
        coeff1 = 1 / (Tc*sqtTr);
        coeff2 = sqtTr - 1;
        k = kappa_vec_[i];
        dalphadT_vec_Curr_[i] = coeff1 *(k* k*coeff2 - k);
    }
    //Calculate mixture derivative
    for (size_t i = 0; i < m_kk; i++) {
        size_t counter1 = i + m_kk * i;
        for (size_t j = 0; j < m_kk; j++) {
            size_t counter2 = j * m_kk + j;
            temp = 0.5 * sqrt((a_vec_Curr_[counter1] * a_vec_Curr_[counter2]) / (alpha_vec_Curr_[i] * alpha_vec_Curr_[j]));
            daAlphadT += moleFractions_[i] * moleFractions_[j] * temp
                        * (dalphadT_vec_Curr_[j] * alpha_vec_Curr_[i] + dalphadT_vec_Curr_[i] * alpha_vec_Curr_[j]);
        }
    }
    return daAlphadT;
}

double PengRobinson::d2aAlpha_dT2() const
{
    double temp, fac1, fac2, alphaij, alphai, alphaj, d2aAlphadT2 = 0.0, num;
    double k;
    for (size_t i = 0; i < m_kk; i++) {
        size_t counter = i + m_kk * i;
        double Tcrit_i = speciesCritTemperature(a_vec_Curr_[counter], b_vec_Curr_[i]);
        double sqt_Tr = sqrt(temperature() / Tcrit_i); //we need species critical temperature
        double coeff1 = 1 / (Tcrit_i*Tcrit_i*sqt_Tr);
        double coeff2 = sqt_Tr - 1;
        //  Calculate first and second derivatives of alpha for individual species
        k = kappa_vec_[i];
        dalphadT_vec_Curr_[i] = coeff1 *(k* k*coeff2 - k);
        d2alphadT2_[i] = (k*k + k) * coeff1 / (2 * sqt_Tr*sqt_Tr);
    }

    //Calculate mixture derivative
    for (size_t i = 0; i < m_kk; i++) {
        size_t counter1 = i + m_kk * i;
        alphai = alpha_vec_Curr_[i];
        for (size_t j = 0; j < m_kk; j++) {
            size_t counter2 = j + m_kk * j;
            alphaj = alpha_vec_Curr_[j];
            alphaij = alphai * alphaj;
            temp = 0.5 * sqrt((a_vec_Curr_[counter1] * a_vec_Curr_[counter2]) / (alphaij));
            num = (dalphadT_vec_Curr_[j] * alphai + dalphadT_vec_Curr_[i] * alphaj);
            fac1 = -(0.5 / alphaij)*num*num;
            fac2 = alphaj * d2alphadT2_[counter1] + alphai *d2alphadT2_[counter2] + 2 * dalphadT_vec_Curr_[i] * dalphadT_vec_Curr_[j];
            d2aAlphadT2 += moleFractions_[i] * moleFractions_[j] * temp *(fac1 + fac2);
        }
    }
    return d2aAlphadT2;
}

void PengRobinson::calcCriticalConditions(double a, double b,
        double& pc, double& tc, double& vc) const
{
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
    tc = a * omega_b / (b * omega_a * GasConstant);
    pc = omega_b * GasConstant * tc / b;
    vc = omega_vc * GasConstant * tc / pc;
}

int PengRobinson::NicholsCall(double T, double pres, double a, double b, double aAlpha, double Vroot[3]) const
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

    int nSolnValues = NicholsSolve(T, pres, a, b, aAlpha, Vroot, an, bn, cn, dn, tc);

    return nSolnValues;
}

}
