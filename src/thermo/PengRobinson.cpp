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

<<<<<<< HEAD
=======
PengRobinson::PengRobinson(XML_Node& phaseRefRoot, const std::string& id_) :
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

>>>>>>> Removing MFTP suffix from PengRobinson class
void PengRobinson::calculateAlpha(const std::string& species, double a, double b, double w)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
<<<<<<< HEAD
        throw CanteraError("PengRobinson::calculateAlpha",
=======
        throw CanteraError("PengRobinson::setSpeciesCoeffs",
>>>>>>> Removing MFTP suffix from PengRobinson class
            "Unknown species '{}'.", species);
    }

    // Calculate value of kappa (independent of temperature)
    // w is an acentric factor of species and must be specified in the CTI file
<<<<<<< HEAD
    if (isnan(w)){
        throw CanteraError("PengRobinson::calculateAlpha",
            "No acentric factor loaded.");
    } else if (w <= 0.491) {
=======

    if (w <= 0.491) {
>>>>>>> Removing MFTP suffix from PengRobinson class
        kappa_vec_[k] = 0.37464 + 1.54226*w - 0.26992*w*w;
    } else {
        kappa_vec_[k] = 0.374642 + 1.487503*w - 0.164423*w*w + 0.016666*w*w*w;
    }

    //Calculate alpha (temperature dependent interaction parameter)
<<<<<<< HEAD
    double critTemp = speciesCritTemperature(a, b); // critical temperature of individual species
    double sqt_T_r = sqrt(temperature() / critTemp);
=======
    double criTemp = speciesCritTemperature(a, b); // critical temperature of individual species
    double sqt_T_r = sqrt(temperature() / criTemp);
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
double PengRobinson::cp_mole() const
{
    _updateReferenceStateThermo();
    double T = temperature();
=======
double PengRobinson::enthalpy_mole() const
{
    _updateReferenceStateThermo();
    double h_ideal = RT() * mean_X(m_h0_RT);
    double h_nonideal = hresid();
    return h_ideal + h_nonideal;
}

double PengRobinson::entropy_mole() const
{
    _updateReferenceStateThermo();
    double sr_ideal = GasConstant * (mean_X(m_s0_R) - sum_xlogx()
                                                - std::log(pressure()/refPressure()));
    double sr_nonideal = sresid();
    return sr_ideal + sr_nonideal;
}

double PengRobinson::cp_mole() const
{
    _updateReferenceStateThermo();
    double TKelvin = temperature();
>>>>>>> Removing MFTP suffix from PengRobinson class
    double mv = molarVolume();
    double vpb = mv + (1 + M_SQRT2)*m_b_current;
    double vmb = mv + (1 - M_SQRT2)*m_b_current;
    pressureDerivatives();
    double cpref = GasConstant * mean_X(m_cp0_R);
    double dHdT_V = cpref + mv * dpdT_ - GasConstant
<<<<<<< HEAD
                        + 1.0 / (2.0 * M_SQRT2 *m_b_current) * log(vpb / vmb) * T *d2aAlpha_dT2();
    return dHdT_V - (mv + T * dpdT_ / dpdV_) * dpdT_;
=======
                        + 1.0 / (2.0 * M_SQRT2 *m_b_current) * log(vpb / vmb) * TKelvin *d2aAlpha_dT2();
    return dHdT_V - (mv + TKelvin * dpdT_ / dpdV_) * dpdT_;
>>>>>>> Removing MFTP suffix from PengRobinson class
}

double PengRobinson::cv_mole() const
{
    _updateReferenceStateThermo();
<<<<<<< HEAD
    double T = temperature();
    pressureDerivatives();
    return (cp_mole() + T* dpdT_* dpdT_ / dpdV_);
=======
    double TKelvin = temperature();
    pressureDerivatives();
    return (cp_mole() + TKelvin* dpdT_* dpdT_ / dpdV_);
>>>>>>> Removing MFTP suffix from PengRobinson class
}

double PengRobinson::pressure() const
{
    _updateReferenceStateThermo();
    //  Get a copy of the private variables stored in the State object
<<<<<<< HEAD
    double T = temperature();
    double mv = molarVolume();
    double denom = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    double pp = GasConstant * T / (mv - m_b_current) - m_aAlpha_current / denom;
    return pp;
}

=======
    double TKelvin = temperature();
    double mv = molarVolume();
    double denom = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    double pp = GasConstant * TKelvin / (mv - m_b_current) - m_aAlpha_current / denom;
    return pp;
}

void PengRobinson::calcDensity()
{
    // Calculate the molarVolume of the solution (m**3 kmol-1)
    const double* const dtmp = moleFractdivMMW();
    getPartialMolarVolumes(m_tmpV.data());
    double invDens = dot(m_tmpV.begin(), m_tmpV.end(), dtmp);

    // Set the density in the parent State object directly, by calling the
    // Phase::setDensity() function.
    Phase::setDensity(1.0/invDens);
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

void PengRobinson::getActivityConcentrations(double* c) const
{
    getActivityCoefficients(c);
    double p_RT = pressure() / RT();
    for (size_t k = 0; k < m_kk; k++) {
        c[k] *= moleFraction(k)* p_RT;
    }
}

>>>>>>> Removing MFTP suffix from PengRobinson class
double PengRobinson::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void PengRobinson::getActivityCoefficients(double* ac) const
{
    double mv = molarVolume();
<<<<<<< HEAD
    //double T = temperature();
=======
    double T = temperature();
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
    double T = temperature();
=======
    double TKelvin = temperature();
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
    double fac = T * daAlphadT - m_aAlpha_current;

    pressureDerivatives();
    double fac2 = mv + T * dpdT_ / dpdV_;
=======
    double fac = TKelvin * daAlphadT - m_aAlpha_current;

    pressureDerivatives();
    double fac2 = mv + TKelvin * dpdT_ / dpdV_;
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
        sbar[k] += GasConstant * log(GasConstant * T / (refP * mv))
=======
        sbar[k] += GasConstant * log(GasConstant * TKelvin / (refP * mv))
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
        vector_fp spCoeff{ NAN, NAN, NAN };
=======
        vector_fp spCoeff{ NAN, NAN };
>>>>>>> Removing MFTP suffix from PengRobinson class

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
<<<<<<< HEAD
                double T_crit = 0.0, P_crit = 0.0, w_ac = 0.0;
=======
                double T_crit, P_crit;
>>>>>>> Removing MFTP suffix from PengRobinson class

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
<<<<<<< HEAD
                if (acNodeDoc.hasChild("omega")) {
                    vParams = 0.0;
                    XML_Node& xmlChildCoeff = acNodeDoc.child("omega");
                    if (xmlChildCoeff.hasChild("value")) {
                        std::string acentric_factor = xmlChildCoeff.attrib("value");
                        vParams = strSItoDbl(acentric_factor);
                    }
                    w_ac = vParams;

                }
=======
>>>>>>> Removing MFTP suffix from PengRobinson class

                //Assuming no temperature dependence
                spCoeff[0] = omega_a * (GasConstant* GasConstant) * (T_crit* T_crit) / P_crit; //coeff a
                spCoeff[1] = omega_b * GasConstant * T_crit / P_crit; // coeff b
<<<<<<< HEAD
                spCoeff[2] = w_ac; // acentric factor
=======
>>>>>>> Removing MFTP suffix from PengRobinson class
                break;
            }
        }
        return spCoeff;
}

<<<<<<< HEAD
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
            if (a_coeff_vec(0, k + m_kk * k) == 0.0) {
                // coeffs[0] = a0, coeffs[1] = b;
                vector<double> coeffs = getCoeff(item.first);

                // Check if species was found in the database of critical
                // properties, and assign the results
                if (!isnan(coeffs[0])) {
                    // Assuming no temperature dependence (i.e. a1 = 0)
                    calculateAlpha(item.first, coeffs[0], coeffs[1], coeffs[2]);
                    setSpeciesCoeffs(item.first, coeffs[0], coeffs[1], coeffs[2]);
                }
            }
        }
    }
=======
void PengRobinson::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
        if (phaseNode.hasChild("thermo")) {
            XML_Node& thermoNode = phaseNode.child("thermo");
            std::string model = thermoNode["model"];
            if (model != "PengRobinson" && model != "PengRobinson") {
                throw CanteraError("PengRobinson::initThermoXML",
                    "Unknown thermo model : " + model);
            }

            // Go get all of the coefficients and factors in the
            // activityCoefficients XML block
            if (thermoNode.hasChild("activityCoefficients")) {
                XML_Node& acNode = thermoNode.child("activityCoefficients");

                // Count the number of species with parameters provided in the
                //    input file:
                size_t nParams = 0;

                // Loop through the children and read out fluid parameters.  Process
                //   all the pureFluidParameters, first:
                for (size_t i = 0; i < acNode.nChildren(); i++) {
                    XML_Node& xmlACChild = acNode.child(i);
                    if (caseInsensitiveEquals(xmlACChild.name(), "purefluidparameters")) {
                        readXMLPureFluid(xmlACChild);
                        nParams += 1;
                    }
                }

                // If any species exist which have undefined pureFluidParameters,
                // search the database in 'critProperties.xml' to find critical
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
                        // coeffArray[0] = a0, coeffArray[1] = b, coeffArray[2] = w;
                        coeffArray = getCoeff(iName);

                        // Check if species was found in the database of critical properties,
                        // and assign the results
                        if (!isnan(coeffArray[0])) {
                                //Assuming no temperature dependence (i,e a1 = 0)
                                setSpeciesCoeffs(iName, coeffArray[0], 0.0, coeffArray[1]);
                        }
                    }
                }

                // Loop back through the "activityCoefficients" children and process the
                // crossFluidParameters in the XML tree:
                for (size_t i = 0; i < acNode.nChildren(); i++) {
                    XML_Node& xmlACChild = acNode.child(i);
                    if (caseInsensitiveEquals(xmlACChild.name(), "crossfluidparameters")) {
                        readXMLCrossFluid(xmlACChild);
                    }
                }
            }
        }

        MixtureFugacityTP::initThermoXML(phaseNode, id);
}

void PengRobinson::readXMLPureFluid(XML_Node& pureFluidParam)
{
        string xname = pureFluidParam.name();
        if (xname != "pureFluidParameters") {
            throw CanteraError("PengRobinson::readXMLPureFluid",
                "Incorrect name for processing this routine: " + xname);
        }

        double a0 = 0.0;
        double a1 = 0.0;
        double b = 0.0;
        double w = 0.0;
        for (size_t iChild = 0; iChild < pureFluidParam.nChildren(); iChild++) {
            XML_Node& xmlChild = pureFluidParam.child(iChild);
            string nodeName = toLowerCopy(xmlChild.name());

            if (nodeName == "a_coeff") {
                vector_fp vParams;
                string iModel = toLowerCopy(xmlChild.attrib("model"));
                getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");

                if (vParams.size() == 1) {
                    a0 = vParams[0];
                } else if (vParams.size() == 2) {
                    a0 = vParams[0];
                    a1 = vParams[1];
                } else {
                    throw CanteraError("PengRobinson::readXMLPureFluid",
                        "unknown model or incorrect number of parameters");
                }
            } else if (nodeName == "b_coeff") {
                b = getFloatCurrent(xmlChild, "toSI");
            } else if (nodeName == "acentric_factor") {
                w = getFloatCurrent(xmlChild);
            }
        }
        calculateAlpha(pureFluidParam.attrib("species"), a0, b, w);
        setSpeciesCoeffs(pureFluidParam.attrib("species"), a0, b, w);
}

void PengRobinson::readXMLCrossFluid(XML_Node& CrossFluidParam)
{
        string xname = CrossFluidParam.name();
        if (xname != "crossFluidParameters") {
            throw CanteraError("PengRobinson::readXMLCrossFluid",
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
                    throw CanteraError("PengRobinson::readXMLCrossFluid",
                        "unknown model ({}) or wrong number of parameters ({})",
                        iModel, vParams.size());
                }
            }
        }
}

void PengRobinson::setParametersFromXML(const XML_Node& thermoNode)
{
    MixtureFugacityTP::setParametersFromXML(thermoNode);
>>>>>>> Removing MFTP suffix from PengRobinson class
}

double PengRobinson::sresid() const
{
    double molarV = molarVolume();
    double hh = m_b_current / molarV;
    double zz = z();
    double alpha_1 = daAlpha_dT();
<<<<<<< HEAD
=======
    double T = temperature();
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
double PengRobinson::liquidVolEst(double T, double& presGuess) const
=======
double PengRobinson::liquidVolEst(double TKelvin, double& presGuess) const
>>>>>>> Removing MFTP suffix from PengRobinson class
{
    double v = m_b_current * 1.1;
    double atmp;
    double btmp;
    double aAlphatmp;
<<<<<<< HEAD
    calculateAB(T, atmp, btmp, aAlphatmp);
    double pres = std::max(psatEst(T), presGuess);
=======
    calculateAB(TKelvin, atmp, btmp, aAlphatmp);
    double pres = std::max(psatEst(TKelvin), presGuess);
>>>>>>> Removing MFTP suffix from PengRobinson class
    double Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
<<<<<<< HEAD
                int nsol = NicholsCall(T, pres, atmp, btmp, aAlphatmp, Vroot);
=======
                int nsol = NicholsSolve(TKelvin, pres, atmp, btmp, aAlphatmp, Vroot);
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
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
=======
double PengRobinson::densityCalc(double TKelvin, double presPa, int phaseRequested, double rhoGuess)
{
    // It's necessary to set the temperature so that m_aAlpha_current is set correctly.
    setTemperature(TKelvin);
    double tcrit = critTemperature();
    double mmw = meanMolecularWeight();
    if (rhoGuess == -1.0) {
        if (phaseRequested != FLUID_GAS) {
            if (TKelvin > tcrit) {
                rhoGuess = presPa * mmw / (GasConstant * TKelvin);
            } else {
                if (phaseRequested == FLUID_UNDEFINED || phaseRequested == FLUID_SUPERCRIT) {
                    rhoGuess = presPa * mmw / (GasConstant * TKelvin);
                } else if (phaseRequested >= FLUID_LIQUID_0) {
                    double lqvol = liquidVolEst(TKelvin, presPa);
                    rhoGuess = mmw / lqvol;
                }
            }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to the routine
            rhoGuess = presPa * mmw / (GasConstant * TKelvin);
        }
    }

    double volGuess = mmw / rhoGuess;
    NSolns_ = NicholsSolve(TKelvin, presPa, m_a_current, m_b_current, m_aAlpha_current, Vroot_);
>>>>>>> Removing MFTP suffix from PengRobinson class

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
<<<<<<< HEAD
        } else if (T > tcrit) {
=======
        } else if (TKelvin > tcrit) {
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
    int nsol = NicholsCall(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
=======
    int nsol = NicholsSolve(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
    int nsol = NicholsCall(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
=======
    int nsol = NicholsSolve(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
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
=======
double PengRobinson::pressureCalc(double TKelvin, double molarVol) const
{
    double den = molarVol * molarVol + 2 * molarVol * m_b_current - m_b_current * m_b_current;
    double pres = GasConstant * TKelvin / (molarVol - m_b_current) - m_aAlpha_current / den;
    return pres;
}

double PengRobinson::dpdVCalc(double TKelvin, double molarVol, double& presCalc) const
{
    double den = molarVol * molarVol + 2 * molarVol * m_b_current - m_b_current * m_b_current;
    presCalc = GasConstant * TKelvin / (molarVol - m_b_current) - m_aAlpha_current/ den;

    double vpb = molarVol + m_b_current;
    double vmb = molarVol - m_b_current;
    double dpdv = -GasConstant * TKelvin / (vmb * vmb) + 2 *m_aAlpha_current * vpb / (den*den);
>>>>>>> Removing MFTP suffix from PengRobinson class
    return dpdv;
}

void PengRobinson::pressureDerivatives() const
{
<<<<<<< HEAD
    double T = temperature();
    double mv = molarVolume();
    double pres;

    dpdV_ = dpdVCalc(T, mv, pres);
=======
    double TKelvin = temperature();
    double mv = molarVolume();
    double pres;

    dpdV_ = dpdVCalc(TKelvin, mv, pres);
>>>>>>> Removing MFTP suffix from PengRobinson class
    double vmb = mv - m_b_current;
    double den = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
    dpdT_ = (GasConstant / vmb - daAlpha_dT() / den);
}

void PengRobinson::updateMixingExpressions()
{
<<<<<<< HEAD
    double temp = temperature();
    //Update aAlpha_i
    double sqt_alpha;

    // Update indiviual alpha
    for (size_t j = 0; j < m_kk; j++) {
        size_t counter = j * m_kk + j;
        double critTemp_j = speciesCritTemperature(a_vec_Curr_[counter],b_vec_Curr_[j]);
        sqt_alpha = 1 + kappa_vec_[j] * (1 - sqrt(temp / critTemp_j));
=======
    updateAB();
}

void PengRobinson::updateAB()
{
    double temp = temperature();
    //Update aAlpha_i
    double sqt_alpha;
    double criTemp = critTemperature();
    double sqt_T_reduced = sqrt(temp / criTemp);

    // Update indiviual alpha
    for (size_t j = 0; j < m_kk; j++) {
        sqt_alpha = 1 + kappa_vec_[j] * (1 - sqt_T_reduced);
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
            m_aAlpha_current += aAlpha_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
=======
                        m_aAlpha_current += aAlpha_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
    double temp, fac1, fac2, alphaij, alphai, alphaj, d2aAlphadT2 = 0.0, num;
    double k;
    for (size_t i = 0; i < m_kk; i++) {
        size_t counter = i + m_kk * i;
        double Tcrit_i = speciesCritTemperature(a_vec_Curr_[counter], b_vec_Curr_[i]);
        double sqt_Tr = sqrt(temperature() / Tcrit_i); //we need species critical temperature
        double coeff1 = 1 / (Tcrit_i*Tcrit_i*sqt_Tr);
        double coeff2 = sqt_Tr - 1;
        //  Calculate first and second derivatives of alpha for individual species
=======
    double daAlphadT = 0.0, temp, fac1, fac2, alphaij, alphai, alphaj, d2aAlphadT2 = 0.0, num;
    double k;
    double sqt_Tr = sqrt(temperature() / critTemperature()); //we need species critical temperature
    double coeff1 = 1 / (critTemperature()*critTemperature()*sqt_Tr);
    double coeff2 = sqt_Tr - 1;
    for (size_t i = 0; i < m_kk; i++) {
        //  Calculate first and second derivatives of alpha for individual species
        size_t counter = i + m_kk * i;
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
int PengRobinson::NicholsCall(double T, double pres, double a, double b, double aAlpha, double Vroot[3]) const
{
    // Derive the coefficients of the cubic polynomial (in terms of molar volume v) to solve.
    double bsqr = b * b;
    double RT_p = GasConstant * T / pres;
=======
int PengRobinson::NicholsSolve(double TKelvin, double pres, double a, double b, double aAlpha,
                                   double Vroot[3]) const
{
    double tmp;
    fill_n(Vroot, 3, 0.0);
    if (TKelvin <= 0.0) {
        throw CanteraError("PengRobinson::NicholsSolve()", "negative temperature T = {}", TKelvin);
    }

    // Derive the coefficients of the cubic polynomial (in terms of molar volume v) to solve.
    double bsqr = b * b;
    double RT_p = GasConstant * TKelvin / pres;
>>>>>>> Removing MFTP suffix from PengRobinson class
    double aAlpha_p = aAlpha / pres;
    double an = 1.0;
    double bn = (b - RT_p);
    double cn = -(2 * RT_p * b - aAlpha_p + 3 * bsqr);
    double dn = (bsqr * RT_p + bsqr * b - aAlpha_p * b);

    double tc = a * omega_b / (b * omega_a * GasConstant);
<<<<<<< HEAD

    int nSolnValues = NicholsSolve(T, pres, a, b, aAlpha, Vroot, an, bn, cn, dn, tc);

=======
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    // Derive the center of the cubic, x_N
    double xN = - bn /(3 * an);

    // Derive the value of delta**2. This is a key quantity that determines the number of turning points
    double delta2 = (bn * bn - 3 * an * cn) / (9 * an * an); 
    double delta = 0.0;

    // Calculate a couple of ratios
    // Cubic equation in z : z^3 - (1-B) z^2 + (A -2B -3B^2)z - (AB- B^2- B^3) = 0
    double ratio1 = 3.0 * an * cn / (bn * bn);
    double ratio2 = pres * b / (GasConstant * TKelvin); // B
    if (fabs(ratio1) < 1.0E-7) {
        double ratio3 = aAlpha / (GasConstant * TKelvin) * pres / (GasConstant * TKelvin); // A
        if (fabs(ratio2) < 1.0E-5 && fabs(ratio3) < 1.0E-5) {
            // A and B terms in cubic equation for z are almost zero, then z is near to 1
            double zz = 1.0;
            for (int i = 0; i < 10; i++) {
                double znew = zz / (zz - ratio2) - ratio3 / (zz + ratio1);
                double deltaz = znew - zz;
                zz = znew;
                if (fabs(deltaz) < 1.0E-14) {
                    break;
                }
            }
            double v = zz * GasConstant * TKelvin / pres;
            Vroot[0] = v;
            return 1;
        }
    }

    int nSolnValues; // Represents number of solutions to the cubic equation
    double h2 = 4. * an * an * delta2 * delta2 * delta2; // h^2
    if (delta2 > 0.0) {
        delta = sqrt(delta2);
    }

    double h = 2.0 * an * delta * delta2;
    double yN = 2.0 * bn * bn * bn / (27.0 * an * an) - bn * cn / (3.0 * an) + dn; // y_N term
    double disc = yN * yN - h2; // discriminant

    //check if y = h
    if (fabs(fabs(h) - fabs(yN)) < 1.0E-10) {
        if (disc > 1e-10) {
            throw CanteraError("PengRobinson::NicholsSolve()", "value of yN and h are too high, unrealistic roots may be obtained");
        } 
        disc = 0.0;
    }

    if (disc < -1e-14) {
        // disc<0 then we have three distinct roots.
        nSolnValues = 3;
    } else if (fabs(disc) < 1e-14) {
        // disc=0 then we have two distinct roots (third one is repeated root)
        nSolnValues = 2;
        // We are here as p goes to zero.
    } else if (disc > 1e-14) {
        // disc> 0 then we have one real root.
        nSolnValues = 1;
    }

    // One real root -> have to determine whether gas or liquid is the root
    if (disc > 0.0) {
        double tmpD = sqrt(disc);
        double tmp1 = (- yN + tmpD) / (2.0 * an);
        double sgn1 = 1.0;
        if (tmp1 < 0.0) {
            sgn1 = -1.0;
            tmp1 = -tmp1;
        }
        double tmp2 = (- yN - tmpD) / (2.0 * an);
        double sgn2 = 1.0;
        if (tmp2 < 0.0) {
            sgn2 = -1.0;
            tmp2 = -tmp2;
        }
        double p1 = pow(tmp1, 1./3.);
        double p2 = pow(tmp2, 1./3.);
        double alpha = xN + sgn1 * p1 + sgn2 * p2;
        Vroot[0] = alpha;
        Vroot[1] = 0.0;
        Vroot[2] = 0.0;
    } else if (disc < 0.0) {
        // Three real roots alpha, beta, gamma are obtained.
        double val = acos(-yN / h);
        double theta = val / 3.0;
        double twoThirdPi = 2. * Pi / 3.;
        double alpha = xN + 2. * delta * cos(theta);
        double beta = xN + 2. * delta * cos(theta + twoThirdPi);
        double gamma = xN + 2. * delta * cos(theta + 2.0 * twoThirdPi);
        Vroot[0] = beta;
        Vroot[1] = gamma;
        Vroot[2] = alpha;

        for (int i = 0; i < 3; i++) {
            tmp = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(tmp) > 1.0E-4) {
                for (int j = 0; j < 3; j++) {
                    if (j != i && fabs(Vroot[i] - Vroot[j]) < 1.0E-4 * (fabs(Vroot[i]) + fabs(Vroot[j]))) {
                        writelog("PengRobinson::NicholsSolve(T = {}, p = {}):"
                                 " WARNING roots have merged: {}, {}\n",
                                 TKelvin, pres, Vroot[i], Vroot[j]);
                    }
                }
            }
        }
    } else if (disc == 0.0) {
        //Three equal roots are obtained, i.e. alpha = beta = gamma
        if (yN < 1e-18 && h < 1e-18) {
            // yN = 0.0 and h = 0 i.e. disc = 0
            Vroot[0] = xN;
            Vroot[1] = xN;
            Vroot[2] = xN;
        } else {
            // h and yN need to figure out whether delta^3 is positive or negative
            if (yN > 0.0) {
                tmp = pow(yN/(2*an), 1./3.);
                // In this case, tmp and delta must be equal. 
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("PengRobinson::NicholsSolve()", "Inconsistancy in cubic solver : solver is bad conditioned.");
                }
                Vroot[1] = xN + delta;
                Vroot[0] = xN - 2.0*delta; // liquid phase root
            } else {
                tmp = pow(yN/(2*an), 1./3.);
                // In this case, tmp and delta must be equal. 
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("PengRobinson::NicholsSolve()", "Inconsistancy in cubic solver : solver is bad conditioned.");
                }
                delta = -delta;
                Vroot[0] = xN + delta;
                Vroot[1] = xN - 2.0*delta; // gas phase root
            }
        }
    }

    // Find an accurate root, since there might be a heavy amount of roundoff error due to bad conditioning in this solver.
    double res, dresdV = 0.0;
    for (int i = 0; i < nSolnValues; i++) {
        for (int n = 0; n < 20; n++) {
            res = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res) < 1.0E-14) { // accurate root is obtained
                break;
            }
            dresdV = 3.0 * an * Vroot[i] * Vroot[i] + 2.0 * bn * Vroot[i] + cn;     // derivative of the residual
            double del = - res / dresdV;
            Vroot[i] += del;
            if (fabs(del) / (fabs(Vroot[i]) + fabs(del)) < 1.0E-14) {
                break;
            }
            double res2 = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res2) < fabs(res)) {
                continue;
            } else { 
                Vroot[i] -= del;        // Go back to previous value of Vroot.
                Vroot[i] += 0.1 * del;  // under-relax by 0.1
            }
        }
        if ((fabs(res) > 1.0E-14) && (fabs(res) > 1.0E-14 * fabs(dresdV) * fabs(Vroot[i]))) {
            writelog("PengRobinson::NicholsSolve(T = {}, p = {}): "
                "WARNING root didn't converge V = {}", TKelvin, pres, Vroot[i]);
            writelogendl();
        }
    }

    if (nSolnValues == 1) {
        if (TKelvin > tc) {
            if (Vroot[0] < vc) {
                // Liquid phase root
                nSolnValues = -1;
            }
        } else {
            if (Vroot[0] < xN) {
                nSolnValues = -1;
            }
        }
    } else {
        if (nSolnValues == 2 && delta > 1e-14) {
            nSolnValues = -2;
        }
    }
>>>>>>> Removing MFTP suffix from PengRobinson class
    return nSolnValues;
}

}
