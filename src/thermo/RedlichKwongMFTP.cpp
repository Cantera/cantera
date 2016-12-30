//! @file RedlichKwongMFTP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

const doublereal RedlichKwongMFTP::omega_a = 4.27480233540E-01;
const doublereal RedlichKwongMFTP::omega_b = 8.66403499650E-02;
const doublereal RedlichKwongMFTP::omega_vc = 3.33333333333333E-01;

RedlichKwongMFTP::RedlichKwongMFTP() :
    m_standardMixingRules(0),
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    Vroot_{0.0, 0.0, 0.0},
    dpdV_(0.0),
    dpdT_(0.0)
{
}

RedlichKwongMFTP::RedlichKwongMFTP(const std::string& infile, const std::string& id_) :
    m_standardMixingRules(0),
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    Vroot_{0.0, 0.0, 0.0},
    dpdV_(0.0),
    dpdT_(0.0)
{
    initThermoFile(infile, id_);
}

RedlichKwongMFTP::RedlichKwongMFTP(XML_Node& phaseRefRoot, const std::string& id_) :
    m_standardMixingRules(0),
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    Vroot_{0.0, 0.0, 0.0},
    dpdV_(0.0),
    dpdT_(0.0)
{
    importPhase(phaseRefRoot, this);
}

RedlichKwongMFTP::RedlichKwongMFTP(const RedlichKwongMFTP& b) :
    m_standardMixingRules(0),
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    *this = b;
}

RedlichKwongMFTP& RedlichKwongMFTP::operator=(const RedlichKwongMFTP& b)
{
    if (&b != this) {
        // Mostly, this is a passthrough to the underlying assignment operator
        // for the ThermoPhae parent object.
        MixtureFugacityTP::operator=(b);

        // However, we have to handle data that we own.
        m_standardMixingRules = b.m_standardMixingRules;
        m_formTempParam = b.m_formTempParam;
        m_b_current = b.m_b_current;
        m_a_current = b.m_a_current;
        a_vec_Curr_ = b.a_vec_Curr_;
        b_vec_Curr_ = b.b_vec_Curr_;
        a_coeff_vec = b.a_coeff_vec;

        m_pc_Species = b.m_pc_Species;
        m_tc_Species = b.m_tc_Species;
        m_vc_Species = b.m_vc_Species;
        NSolns_ = b.NSolns_;
        Vroot_[0] = b.Vroot_[0];
        Vroot_[1] = b.Vroot_[1];
        Vroot_[2] = b.Vroot_[2];
        m_pp = b.m_pp;
        m_tmpV = b.m_tmpV;
        m_partialMolarVolumes = b.m_partialMolarVolumes;
        dpdV_ = b.dpdV_;
        dpdT_ = b.dpdT_;
        dpdni_ = b.dpdni_;
    }
    return *this;
}

ThermoPhase* RedlichKwongMFTP::duplMyselfAsThermoPhase() const
{
    return new RedlichKwongMFTP(*this);
}

int RedlichKwongMFTP::eosType() const
{
    warn_deprecated("RedlichKwongMFTP::eosType",
                    "To be removed after Cantera 2.3.");
    return cRedlichKwongMFTP;
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
                                          - sum_xlogx() - std::log(pressure()/m_spthermo->refPressure()));
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
    if (fabs(pp -m_Pcurrent) > 1.0E-5 * fabs(m_Pcurrent)) {
        throw CanteraError(" RedlichKwongMFTP::pressure()", "setState broken down, maybe");
    }

    return m_Pcurrent;
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
        doublereal denom = (m_Pcurrent + RT() * m_b_current/(vmb * vmb) - m_a_current / (sqt * vpb * vpb)
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
        b_vec_Curr_.push_back(0.0);

        a_coeff_vec.resize(2, m_kk * m_kk, 0.0);

        m_pc_Species.push_back(0.0);
        m_tc_Species.push_back(0.0);
        m_vc_Species.push_back(0.0);

        m_pp.push_back(0.0);
        m_tmpV.push_back(0.0);
        m_partialMolarVolumes.push_back(0.0);
        dpdni_.push_back(0.0);
    }
    return added;
}

void RedlichKwongMFTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    // Check the model parameter for the Redlich-Kwong equation of state
    // two are allowed
    //       RedlichKwong        mixture of species, each of which are RK fluids
    //       RedlichKwongMFTP    mixture of species with cross term coefficients
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thermoNode = phaseNode.child("thermo");
        std::string model = thermoNode["model"];
        if (model == "RedlichKwong") {
            m_standardMixingRules = 1;
        } else if (model == "RedlichKwongMFTP") {
            m_standardMixingRules = 0;
        } else {
            throw CanteraError("RedlichKwongMFTP::initThermoXML",
                               "Unknown thermo model : " + model);
        }

        // Go get all of the coefficients and factors in the
        // activityCoefficients XML block
        XML_Node* acNodePtr = 0;
        if (thermoNode.hasChild("activityCoefficients")) {
            XML_Node& acNode = thermoNode.child("activityCoefficients");
            acNodePtr = &acNode;
            size_t nC = acNode.nChildren();

            // Loop through the children getting multiple instances of
            // parameters
            for (size_t i = 0; i < nC; i++) {
                XML_Node& xmlACChild = acNodePtr->child(i);
                if (ba::iequals(xmlACChild.name(), "purefluidparameters")) {
                    readXMLPureFluid(xmlACChild);
                }
            }
            if (m_standardMixingRules == 1) {
                applyStandardMixingRules();
            }

            // Loop through the children getting multiple instances of
            // parameters
            for (size_t i = 0; i < nC; i++) {
                XML_Node& xmlACChild = acNodePtr->child(i);
                if (ba::iequals(xmlACChild.name(), "crossfluidparameters")) {
                    readXMLCrossFluid(xmlACChild);
                }
            }
        }
    }

    for (size_t i = 0; i < m_kk; i++) {
        double a0coeff = a_coeff_vec(0, i*m_kk + i);
        double aTcoeff = a_coeff_vec(1, i*m_kk + i);
        double ai = a0coeff + aTcoeff * 500.;
        double bi = b_vec_Curr_[i];
        calcCriticalConditions(ai, bi, a0coeff, aTcoeff, m_pc_Species[i], m_tc_Species[i], m_vc_Species[i]);
    }

    MixtureFugacityTP::initThermoXML(phaseNode, id);
}

void RedlichKwongMFTP::readXMLPureFluid(XML_Node& pureFluidParam)
{
    vector_fp vParams;
    string xname = pureFluidParam.name();
    if (xname != "pureFluidParameters") {
        throw CanteraError("RedlichKwongMFTP::readXMLPureFluid",
                           "Incorrect name for processing this routine: " + xname);
    }

    // Read the species. Find the index of the species in the current phase.
    // It's not an error to not find the species
    string iName = pureFluidParam.attrib("species");
    if (iName == "") {
        throw CanteraError("RedlichKwongMFTP::readXMLPureFluid", "no species attribute");
    }
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    size_t counter = iSpecies + m_kk * iSpecies;
    size_t nParamsExpected, nParamsFound;
    size_t num = pureFluidParam.nChildren();
    for (size_t iChild = 0; iChild < num; iChild++) {
        XML_Node& xmlChild = pureFluidParam.child(iChild);
        string nodeName = ba::to_lower_copy(xmlChild.name());

        if (nodeName == "a_coeff") {
            string iModel = ba::to_lower_copy(xmlChild.attrib("model"));
            if (iModel == "constant") {
                nParamsExpected = 1;
            } else if (iModel == "linear_a") {
                nParamsExpected = 2;
                if (m_formTempParam == 0) {
                    m_formTempParam = 1;
                }
            } else {
                throw CanteraError("RedlichKwongMFTP::readXMLPureFluid", "unknown model");
            }

            getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");
            nParamsFound = vParams.size();
            if (nParamsFound != nParamsExpected) {
                throw CanteraError("RedlichKwongMFTP::readXMLPureFluid(for a_coeff" + iName + ")",
                                   "wrong number of params found");
            }

            for (size_t i = 0; i < nParamsFound; i++) {
                a_coeff_vec(i, counter) = vParams[i];
            }
        } else if (nodeName == "b_coeff") {
            getFloatArray(xmlChild, vParams, true, "m3/kmol", "b_coeff");
            nParamsFound = vParams.size();
            if (nParamsFound != 1) {
                throw CanteraError("RedlichKwongMFTP::readXMLPureFluid(for b_coeff" + iName + ")",
                                   "wrong number of params found");
            }
            b_vec_Curr_[iSpecies] = vParams[0];
        }
    }
}

void RedlichKwongMFTP::applyStandardMixingRules()
{
    int nParam = 2;
    for (size_t i = 0; i < m_kk; i++) {
        size_t icounter = i + m_kk * i;
        for (size_t j = 0; j < m_kk; j++) {
            if (i != j) {
                size_t counter = i + m_kk * j;
                size_t jcounter = j + m_kk * j;
                for (int n = 0; n < nParam; n++) {
                    a_coeff_vec(n, counter) = sqrt(a_coeff_vec(n, icounter) * a_coeff_vec(n, jcounter));
                }
            }
        }
    }
}

void RedlichKwongMFTP::readXMLCrossFluid(XML_Node& CrossFluidParam)
{
    vector_fp vParams;
    string xname = CrossFluidParam.name();
    if (xname != "crossFluidParameters") {
        throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid",
                           "Incorrect name for processing this routine: " + xname);
    }

    // Read the species. Find the index of the species in the current phase.
    // It's not an error to not find the species
    string iName = CrossFluidParam.attrib("species1");
    if (iName == "") {
        throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid", "no species1 attribute");
    }
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    string jName = CrossFluidParam.attrib("species2");
    if (iName == "") {
        throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid", "no species2 attribute");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }

    size_t counter = iSpecies + m_kk * jSpecies;
    size_t counter0 = jSpecies + m_kk * iSpecies;
    size_t nParamsExpected, nParamsFound;
    size_t num = CrossFluidParam.nChildren();
    for (size_t iChild = 0; iChild < num; iChild++) {
        XML_Node& xmlChild = CrossFluidParam.child(iChild);
        string nodeName = ba::to_lower_copy(xmlChild.name());

        if (nodeName == "a_coeff") {
            string iModel = ba::to_lower_copy(xmlChild.attrib("model"));
            if (iModel == "constant") {
                nParamsExpected = 1;
            } else if (iModel == "linear_a") {
                nParamsExpected = 2;
                if (m_formTempParam == 0) {
                    m_formTempParam = 1;
                }
            } else {
                throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid", "unknown model");
            }

            getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");
            nParamsFound = vParams.size();
            if (nParamsFound != nParamsExpected) {
                throw CanteraError("RedlichKwongMFTP::readXMLCrossFluid(for a_coeff" + iName + ")",
                                   "wrong number of params found");
            }

            for (size_t i = 0; i < nParamsFound; i++) {
                a_coeff_vec(i, counter) = vParams[i];
                a_coeff_vec(i, counter0) = vParams[i];
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
        throw CanteraError("RedlichKwongMFTP::NicholsSolve()", "neg temperature");
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
            throw CanteraError("NicholsSolve()", "numerical issues");
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
                        writelog("RedlichKwongMFTP::NicholsSolve(T = {}, p = {}):"
                                 " WARNING roots have merged: {}, {}\n",
                                 TKelvin, pres, Vroot[i], Vroot[j]);
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
                    throw CanteraError("RedlichKwongMFTP::NicholsSolve()", "unexpected");
                }
                Vroot[1] = xN + delta;
                Vroot[0] = xN - 2.0*delta; // liquid phase root
            } else {
                tmp = pow(yN/(2*an), 1./3.);
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("RedlichKwongMFTP::NicholsSolve()", "unexpected");
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
            writelog("RedlichKwongMFTP::NicholsSolve(T = {}, p = {}): "
                "WARNING root didn't converge V = {}", TKelvin, pres, Vroot[i]);
            writelogendl();
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
