//! @file RedlichKwongMFTP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/ThermoFactory.h"
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

RedlichKwongMFTP::RedlichKwongMFTP() :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    fill_n(Vroot_, 3, 0.0);
}

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
    size_t counter1 = ki + m_kk * kj;
    size_t counter2 = kj + m_kk * ki;
    a_coeff_vec(0, counter1) = a_coeff_vec(0, counter2) = a0;
    a_coeff_vec(1, counter1) = a_coeff_vec(1, counter2) = a1;
    a_vec_Curr_[counter1] = a_vec_Curr_[counter2] = a0;
}

// ------------Molar Thermodynamic Properties -------------------------

doublereal RedlichKwongMFTP::enthalpy_mole() const
{
    _updateReferenceStateThermo();
    doublereal h_ideal = RT() * mean_X(m_h0_RT);
    doublereal h_nonideal = hresid();
    return h_ideal + h_nonideal;
}

doublereal RedlichKwongMFTP::entropy_mole() const
{
    _updateReferenceStateThermo();
    doublereal sr_ideal = GasConstant * (mean_X(m_s0_R)
                                          - sum_xlogx() - std::log(pressure()/refPressure()));
    doublereal sr_nonideal = sresid();
    return sr_ideal + sr_nonideal;
}

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

void RedlichKwongMFTP::calcDensity()
{
    // Calculate the molarVolume of the solution (m**3 kmol-1)
    const doublereal* const dtmp = moleFractdivMMW();
    getPartialMolarVolumes(m_tmpV.data());
    double invDens = dot(m_tmpV.begin(), m_tmpV.end(), dtmp);

    // Set the density in the parent State object directly, by calling the
    // Phase::setDensity() function.
    Phase::setDensity(1.0/invDens);
}

void RedlichKwongMFTP::setTemperature(const doublereal temp)
{
    Phase::setTemperature(temp);
    _updateReferenceStateThermo();
    updateAB();
}

void RedlichKwongMFTP::compositionChanged()
{
    MixtureFugacityTP::compositionChanged();
    updateAB();
}

void RedlichKwongMFTP::getActivityConcentrations(doublereal* c) const
{
    getActivityCoefficients(c);
    for (size_t k = 0; k < m_kk; k++) {
        c[k] *= moleFraction(k)*pressure()/RT();
    }
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
    getIntEnergy_RT(ubar);
    scale(ubar, ubar+m_kk, ubar, RT());
}

void RedlichKwongMFTP::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    scale(cpbar, cpbar+m_kk, cpbar, GasConstant);
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

doublereal RedlichKwongMFTP::critTemperature() const
{
    double pc, tc, vc;
    double a0 = 0.0;
    double aT = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j <m_kk; j++) {
            size_t counter = i + m_kk * j;
            a0 += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(0, counter);
            aT += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(1, counter);
        }
    }
    calcCriticalConditions(m_a_current, m_b_current, a0, aT, pc, tc, vc);
    return tc;
}

doublereal RedlichKwongMFTP::critPressure() const
{
    double pc, tc, vc;
    double a0 = 0.0;
    double aT = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j <m_kk; j++) {
            size_t counter = i + m_kk * j;
            a0 += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(0, counter);
            aT += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(1, counter);
        }
    }
    calcCriticalConditions(m_a_current, m_b_current, a0, aT, pc, tc, vc);
    return pc;
}

doublereal RedlichKwongMFTP::critVolume() const
{
    double pc, tc, vc;
    double a0 = 0.0;
    double aT = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j <m_kk; j++) {
            size_t counter = i + m_kk * j;
            a0 += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(0, counter);
            aT += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(1, counter);
        }
    }
    calcCriticalConditions(m_a_current, m_b_current, a0, aT, pc, tc, vc);
    return vc;
}

doublereal RedlichKwongMFTP::critCompressibility() const
{
    double pc, tc, vc;
    double a0 = 0.0;
    double aT = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j <m_kk; j++) {
            size_t counter = i + m_kk * j;
            a0 += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(0, counter);
            aT += moleFractions_[i] * moleFractions_[j] * a_coeff_vec(1, counter);
        }
    }
    calcCriticalConditions(m_a_current, m_b_current, a0, aT, pc, tc, vc);
    return pc*vc/tc/GasConstant;
}

doublereal RedlichKwongMFTP::critDensity() const
{
    double pc, tc, vc;
    double a0 = 0.0;
    double aT = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j <m_kk; j++) {
            size_t counter = i + m_kk * j;
            a0 += moleFractions_[i] * moleFractions_[j] *a_coeff_vec(0, counter);
            aT += moleFractions_[i] * moleFractions_[j] *a_coeff_vec(1, counter);
        }
    }
    calcCriticalConditions(m_a_current, m_b_current, a0, aT, pc, tc, vc);
    double mmw = meanMolecularWeight();
    return mmw / vc;
}

void RedlichKwongMFTP::setToEquilState(const doublereal* mu_RT)
{
    double tmp, tmp2;
    _updateReferenceStateThermo();
    getGibbs_RT_ref(m_tmpV.data());

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    doublereal pres = 0.0;
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
        m_tmpV.push_back(0.0);
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
        // 'critProperties.xml' as part of non-XML initialization, so that
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
                // coeffArray[0] = a0, coeffArray[1] = b;
                coeffArray = getCoeff(iName);

                // Check if species was found in the database of critical properties,
                // and assign the results
                if (!isnan(coeffArray[0])) {
                    //Assuming no temperature dependence (i,e a1 = 0)
                    setSpeciesCoeffs(iName, coeffArray[0], 0.0, coeffArray[1]);
                }
            }
        }
    }

    MixtureFugacityTP::initThermoXML(phaseNode, id);
}

void RedlichKwongMFTP::initThermo()
{
    for (auto& item : m_species) {
        // Read a and b coefficients from species 'input' information (i.e. as
        // specified in a YAML input file)
        if (item.second->input.hasKey("equation-of-state")) {
            auto eos = item.second->input["equation-of-state"].getMapWhere(
                "model", "Redlich-Kwong");
            double a0 = 0, a1 = 0;
            if (eos["a"].isScalar()) {
                a0 = eos.convert("a", "Pa*m^6/kmol^2*K^0.5");
            } else {
                auto avec = eos["a"].asVector<AnyValue>(2);
                a0 = eos.units().convert(avec[0], "Pa*m^6/kmol^2*K^0.5");
                a1 = eos.units().convert(avec[1], "Pa*m^6/kmol^2/K^0.5");
            }
            double b = eos.convert("b", "m^3/kmol");
            setSpeciesCoeffs(item.first, a0, a1, b);
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

vector<double> RedlichKwongMFTP::getCoeff(const std::string& iName)
{
    vector_fp spCoeff{NAN, NAN};

    // Get number of species in the database
    // open xml file critProperties.xml
    XML_Node* doc = get_XML_File("critProperties.xml");
    size_t nDatabase = doc->nChildren();

    // Loop through all species in the database and attempt to match supplied
    // species to each. If present, calculate pureFluidParameters a_k and b_k
    // based on crit  properties T_c and P_c:
    for (size_t isp = 0; isp < nDatabase; isp++) {
        XML_Node& acNodeDoc = doc->child(isp);
        std::string iNameLower = toLowerCopy(iName);
        std::string dbName = toLowerCopy(acNodeDoc.attrib("name"));

        // Attempt to match provided specie iName to current database species
        //  dbName:
        if (iNameLower == dbName) {
            // Read from database and calculate a and b coefficients
            double vParams;
            double T_crit=0.;
            double P_crit=0.;

            if (acNodeDoc.hasChild("Tc")) {
                vParams = 0.0;
                XML_Node& xmlChildCoeff = acNodeDoc.child("Tc");
                if (xmlChildCoeff.hasAttrib("value"))
                {
                    std::string critTemp = xmlChildCoeff.attrib("value");
                    vParams = strSItoDbl(critTemp);
                }
                if (vParams <= 0.0) //Assuming that Pc and Tc are non zero.
                {
                    throw CanteraError("RedlichKwongMFTP::getCoeff",
                                       "Critical Temperature must be positive ");
                }
                T_crit = vParams;
            } else {
                throw CanteraError("RedlichKwongMFTP::getCoeff",
                                   "Critical Temperature not in database ");
            }
            if (acNodeDoc.hasChild("Pc")) {
                vParams = 0.0;
                XML_Node& xmlChildCoeff = acNodeDoc.child("Pc");
                if (xmlChildCoeff.hasAttrib("value"))
                {
                    std::string critPressure = xmlChildCoeff.attrib("value");
                    vParams = strSItoDbl(critPressure);
                }
                if (vParams <= 0.0) //Assuming that Pc and Tc are non zero.
                {
                    throw CanteraError("RedlichKwongMFTP::getCoeff",
                                       "Critical Pressure must be positive ");
                }
                P_crit = vParams;
            } else {
                throw CanteraError("RedlichKwongMFTP::getCoeff",
                                   "Critical Pressure not in database ");
            }

            //Assuming no temperature dependence
            spCoeff[0] = omega_a * pow(GasConstant, 2) * pow(T_crit, 2.5) / P_crit; //coeff a
            spCoeff[1] = omega_b * GasConstant * T_crit / P_crit; // coeff b
            break;
        }
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
        int nsol = NicholsSolve(TKelvin, pres, atmp, btmp, Vroot);
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
    doublereal mmw = meanMolecularWeight();
    if (rhoguess == -1.0) {
        if (phaseRequested != FLUID_GAS) {
            if (TKelvin > tcrit) {
                rhoguess = presPa * mmw / (GasConstant * TKelvin);
            } else {
                if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
                    rhoguess = presPa * mmw / (GasConstant * TKelvin);
                } else if (phaseRequested >= FLUID_LIQUID_0) {
                    double lqvol = liquidVolEst(TKelvin, presPa);
                    rhoguess = mmw / lqvol;
                }
            }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to
            // the routine
            rhoguess = presPa * mmw / (GasConstant * TKelvin);
        }
    }

    doublereal volguess = mmw / rhoguess;
    NSolns_ = NicholsSolve(TKelvin, presPa, m_a_current, m_b_current, Vroot_);

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
    int nsol = NicholsSolve(T, pressure(), m_a_current, m_b_current, Vroot);
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
    int nsol = NicholsSolve(T, pressure(), m_a_current, m_b_current, Vroot);
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

doublereal RedlichKwongMFTP::pressureCalc(doublereal TKelvin, doublereal molarVol) const
{
    doublereal sqt = sqrt(TKelvin);
    double pres = GasConstant * TKelvin / (molarVol - m_b_current)
                  - m_a_current / (sqt * molarVol * (molarVol + m_b_current));
    return pres;
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
    updateAB();
}

void RedlichKwongMFTP::updateAB()
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
                    format_to(b, ", {}", speciesName(k));
                } else {
                    format_to(b, "{}", speciesName(k));
                }
            }
        }
        throw CanteraError("RedlichKwongMFTP::updateAB",
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

void RedlichKwongMFTP::calcCriticalConditions(doublereal a, doublereal b, doublereal a0_coeff, doublereal aT_coeff,
        doublereal& pc, doublereal& tc, doublereal& vc) const
{
    if (m_formTempParam != 0) {
        a = a0_coeff;
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
            f = omega_a * b * GasConstant * tc * sqrttc / omega_b - aT_coeff * tc - a0_coeff;
            dfdt = 1.5 * omega_a * b * GasConstant * sqrttc / omega_b - aT_coeff;
            deltatc = - f / dfdt;
            tc += deltatc;
        }
        if (deltatc > 0.1) {
            throw CanteraError("RedlichKwongMFTP::calcCriticalConditions", "didn't converge");
        }
    }

    pc = omega_b * GasConstant * tc / b;
    vc = omega_vc * GasConstant * tc / pc;
}

int RedlichKwongMFTP::NicholsSolve(double TKelvin, double pres, doublereal a, doublereal b,
                                   doublereal Vroot[3]) const
{
    Vroot[0] = 0.0;
    Vroot[1] = 0.0;
    Vroot[2] = 0.0;
    if (TKelvin <= 0.0) {
        throw CanteraError("RedlichKwongMFTP::NicholsSolve", "neg temperature");
    }

    // Derive the coefficients of the cubic polynomial to solve.
    doublereal an = 1.0;
    doublereal bn = - GasConstant * TKelvin / pres;
    doublereal sqt = sqrt(TKelvin);
    doublereal cn = - (GasConstant * TKelvin * b / pres - a/(pres * sqt) + b * b);
    doublereal dn = - (a * b / (pres * sqt));

    double tmp = a * omega_b / (b * omega_a * GasConstant);
    double pp = 2./3.;
    double tc = pow(tmp, pp);
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;
    // Derive the center of the cubic, x_N
    doublereal xN = - bn /(3 * an);

    // Derive the value of delta**2. This is a key quantity that determines the
    // number of turning points
    doublereal delta2 = (bn * bn - 3 * an * cn) / (9 * an * an);
    doublereal delta = 0.0;

    // Calculate a couple of ratios
    doublereal ratio1 = 3.0 * an * cn / (bn * bn);
    doublereal ratio2 = pres * b / (GasConstant * TKelvin);
    if (fabs(ratio1) < 1.0E-7) {
        doublereal ratio3 = a / (GasConstant * sqt) * pres / (GasConstant * TKelvin);
        if (fabs(ratio2) < 1.0E-5 && fabs(ratio3) < 1.0E-5) {
            doublereal zz = 1.0;
            for (int i = 0; i < 10; i++) {
                doublereal znew = zz / (zz - ratio2) - ratio3 / (zz + ratio1);
                doublereal deltaz = znew - zz;
                zz = znew;
                if (fabs(deltaz) < 1.0E-14) {
                    break;
                }
            }
            doublereal v = zz * GasConstant * TKelvin / pres;
            Vroot[0] = v;
            return 1;
        }
    }

    int nSolnValues;
    double h2 = 4. * an * an * delta2 * delta2 * delta2;
    if (delta2 > 0.0) {
        delta = sqrt(delta2);
    }

    doublereal h = 2.0 * an * delta * delta2;
    doublereal yN = 2.0 * bn * bn * bn / (27.0 * an * an) - bn * cn / (3.0 * an) + dn;
    doublereal desc = yN * yN - h2;

    if (fabs(fabs(h) - fabs(yN)) < 1.0E-10) {
        if (desc != 0.0) {
            // this is for getting to other cases
            throw CanteraError("RedlichKwongMFTP::NicholsSolve", "numerical issues");
        }
        desc = 0.0;
    }

    if (desc < 0.0) {
        nSolnValues = 3;
    } else if (desc == 0.0) {
        nSolnValues = 2;
        // We are here as p goes to zero.
    } else if (desc > 0.0) {
        nSolnValues = 1;
    }

    // One real root -> have to determine whether gas or liquid is the root
    if (desc > 0.0) {
        doublereal tmpD = sqrt(desc);
        doublereal tmp1 = (- yN + tmpD) / (2.0 * an);
        doublereal sgn1 = 1.0;
        if (tmp1 < 0.0) {
            sgn1 = -1.0;
            tmp1 = -tmp1;
        }
        doublereal tmp2 = (- yN - tmpD) / (2.0 * an);
        doublereal sgn2 = 1.0;
        if (tmp2 < 0.0) {
            sgn2 = -1.0;
            tmp2 = -tmp2;
        }
        doublereal p1 = pow(tmp1, 1./3.);
        doublereal p2 = pow(tmp2, 1./3.);
        doublereal alpha = xN + sgn1 * p1 + sgn2 * p2;
        Vroot[0] = alpha;
        Vroot[1] = 0.0;
        Vroot[2] = 0.0;
        tmp = an * Vroot[0] * Vroot[0] * Vroot[0] + bn * Vroot[0] * Vroot[0] + cn * Vroot[0] + dn;
    } else if (desc < 0.0) {
        doublereal tmp = - yN/h;
        doublereal val = acos(tmp);
        doublereal theta = val / 3.0;
        doublereal oo = 2. * Pi / 3.;
        doublereal alpha = xN + 2. * delta * cos(theta);
        doublereal beta = xN + 2. * delta * cos(theta + oo);
        doublereal gamma = xN + 2. * delta * cos(theta + 2.0 * oo);
        Vroot[0] = beta;
        Vroot[1] = gamma;
        Vroot[2] = alpha;

        for (int i = 0; i < 3; i++) {
            tmp = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(tmp) > 1.0E-4) {
                for (int j = 0; j < 3; j++) {
                    if (j != i && fabs(Vroot[i] - Vroot[j]) < 1.0E-4 * (fabs(Vroot[i]) + fabs(Vroot[j]))) {
                        warn_user("RedlichKwongMFTP::NicholsSolve",
                            "roots have merged: {}, {} (T = {}, p = {})",
                            Vroot[i], Vroot[j], TKelvin, pres);
                    }
                }
            }
        }
    } else if (desc == 0.0) {
        if (yN == 0.0 && h == 0.0) {
            Vroot[0] = xN;
            Vroot[1] = xN;
            Vroot[2] = xN;
        } else {
            // need to figure out whether delta is pos or neg
            if (yN > 0.0) {
                tmp = pow(yN/(2*an), 1./3.);
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("RedlichKwongMFTP::NicholsSolve", "unexpected");
                }
                Vroot[1] = xN + delta;
                Vroot[0] = xN - 2.0*delta; // liquid phase root
            } else {
                tmp = pow(yN/(2*an), 1./3.);
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("RedlichKwongMFTP::NicholsSolve", "unexpected");
                }
                delta = -delta;
                Vroot[0] = xN + delta;
                Vroot[1] = xN - 2.0*delta; // gas phase root
            }
        }
        for (int i = 0; i < 2; i++) {
            tmp = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
        }
    }

    // Unfortunately, there is a heavy amount of roundoff error due to bad
    // conditioning in this
    double res, dresdV = 0.0;
    for (int i = 0; i < nSolnValues; i++) {
        for (int n = 0; n < 20; n++) {
            res = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res) < 1.0E-14) {
                break;
            }
            dresdV = 3.0 * an * Vroot[i] * Vroot[i] + 2.0 * bn * Vroot[i] + cn;
            double del = - res / dresdV;
            Vroot[i] += del;
            if (fabs(del) / (fabs(Vroot[i]) + fabs(del)) < 1.0E-14) {
                break;
            }
            double res2 = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res2) < fabs(res)) {
                continue;
            } else {
                Vroot[i] -= del;
                Vroot[i] += 0.1 * del;
            }
        }
        if ((fabs(res) > 1.0E-14) && (fabs(res) > 1.0E-14 * fabs(dresdV) * fabs(Vroot[i]))) {
            warn_user("RedlichKwongMFTP::NicholsSolve",
                "root did not converge: V = {} (T = {}, p = {})",
                Vroot[i], TKelvin, pres);
        }
    }

    if (nSolnValues == 1) {
        if (TKelvin > tc) {
            if (Vroot[0] < vc) {
                nSolnValues = -1;
            }
        } else {
            if (Vroot[0] < xN) {
                nSolnValues = -1;
            }
        }
    } else {
        if (nSolnValues == 2 && delta > 0.0) {
            nSolnValues = -2;
        }
    }
    return nSolnValues;
}

}
