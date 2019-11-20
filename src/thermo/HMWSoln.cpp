/**
 *  @file HMWSoln.cpp
 *    Definitions for the HMWSoln ThermoPhase object, which
 *    models concentrated electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::HMWSoln HMWSoln \endlink) .
 *
 * Class HMWSoln represents a concentrated liquid electrolyte phase which obeys
 * the Pitzer formulation for nonideality using molality-based standard states.
 *
 * This version of the code was modified to have the binary Beta2 Pitzer
 * parameter consistent with the temperature expansions used for Beta0,
 * Beta1, and Cphi.(CFJC, SNL)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/electrolytes.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

HMWSoln::HMWSoln() :
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_molalitiesAreCropped(false),
    IMS_X_o_cutoff_(0.2),
    IMS_cCut_(0.05),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0)
{
}

HMWSoln::~HMWSoln()
{
}

HMWSoln::HMWSoln(const std::string& inputFile, const std::string& id_) :
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_molalitiesAreCropped(false),
    IMS_X_o_cutoff_(0.2),
    IMS_cCut_(0.05),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0)
{
    initThermoFile(inputFile, id_);
}

HMWSoln::HMWSoln(XML_Node& phaseRoot, const std::string& id_) :
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_molalitiesAreCropped(false),
    IMS_X_o_cutoff_(0.2),
    IMS_cCut_(0.05),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0)
{
    importPhase(phaseRoot, this);
}

// -------- Molar Thermodynamic Properties of the Solution ---------------

doublereal HMWSoln::enthalpy_mole() const
{
    getPartialMolarEnthalpies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal HMWSoln::relative_enthalpy() const
{
    getPartialMolarEnthalpies(m_tmpV.data());
    double hbar = mean_X(m_tmpV);
    getEnthalpy_RT(m_gamma_tmp.data());
    for (size_t k = 0; k < m_kk; k++) {
        m_gamma_tmp[k] *= RT();
    }
    double h0bar = mean_X(m_gamma_tmp);
    return hbar - h0bar;
}

doublereal HMWSoln::relative_molal_enthalpy() const
{
    double L = relative_enthalpy();
    getMoleFractions(m_tmpV.data());
    double xanion = 0.0;
    size_t kcation = npos;
    double xcation = 0.0;
    size_t kanion = npos;
    for (size_t k = 0; k < m_kk; k++) {
        if (charge(k) > 0.0) {
            if (m_tmpV[k] > xanion) {
                xanion = m_tmpV[k];
                kanion = k;
            }
        } else if (charge(k) < 0.0) {
            if (m_tmpV[k] > xcation) {
                xcation = m_tmpV[k];
                kcation = k;
            }
        }
    }
    if (kcation == npos || kanion == npos) {
        return L;
    }
    double xuse = xcation;
    double factor = 1;
    if (xanion < xcation) {
        xuse = xanion;
        if (charge(kcation) != 1.0) {
            factor = charge(kcation);
        }
    } else {
        if (charge(kanion) != 1.0) {
            factor = charge(kanion);
        }
    }
    xuse = xuse / factor;
    return L / xuse;
}

doublereal HMWSoln::entropy_mole() const
{
    getPartialMolarEntropies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal HMWSoln::gibbs_mole() const
{
    getChemPotentials(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal HMWSoln::cp_mole() const
{
    getPartialMolarCp(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal HMWSoln::cv_mole() const
{
    double kappa_t = isothermalCompressibility();
    double beta = thermalExpansionCoeff();
    double cp = cp_mole();
    double tt = temperature();
    double molarV = molarVolume();
    return cp - beta * beta * tt * molarV / kappa_t;
}

// ------- Mechanical Equation of State Properties ------------------------

void HMWSoln::calcDensity()
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if(cached.validate(temperature(), pressure(), stateMFNumber())) {
        return;
    }

    // Calculate all of the other standard volumes. Note these are constant for
    // now
    getPartialMolarVolumes(m_tmpV.data());
    double dd = meanMolecularWeight() / mean_X(m_tmpV);
    Phase::assignDensity(dd);
}

void HMWSoln::setDensity(const doublereal rho)
{
    warn_deprecated("HMWSoln::setDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setDensity instead");
    double dens_old = density();
    if (rho != dens_old) {
        throw CanteraError("HMWSoln::setDensity",
                           "Density is not an independent variable");
    }
}

void HMWSoln::setMolarDensity(const doublereal rho)
{
    warn_deprecated("HMWSoln::setMolarDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setMolarDensity instead");
    throw CanteraError("HMWSoln::setMolarDensity",
                       "Density is not an independent variable");
}

// ------- Activities and Activity Concentrations

void HMWSoln::getActivityConcentrations(doublereal* c) const
{
    double cs_solvent = standardConcentration();
    getActivities(c);
    c[0] *= cs_solvent;
    if (m_kk > 1) {
        double cs_solute = standardConcentration(1);
        for (size_t k = 1; k < m_kk; k++) {
            c[k] *= cs_solute;
        }
    }
}

doublereal HMWSoln::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    double mvSolvent = m_tmpV[0];
    if (k > 0) {
        return m_Mnaught / mvSolvent;
    }
    return 1.0 / mvSolvent;
}

void HMWSoln::getActivities(doublereal* ac) const
{
    updateStandardStateThermo();

    // Update the molality array, m_molalities(). This requires an update due to
    //   mole fractions
    s_update_lnMolalityActCoeff();

    // Now calculate the array of activities.
    for (size_t k = 1; k < m_kk; k++) {
        ac[k] = m_molalities[k] * exp(m_lnActCoeffMolal_Scaled[k]);
    }
    double xmolSolvent = moleFraction(0);
    ac[0] = exp(m_lnActCoeffMolal_Scaled[0]) * xmolSolvent;
}

void HMWSoln::getUnscaledMolalityActivityCoefficients(doublereal* acMolality) const
{
    updateStandardStateThermo();
    A_Debye_TP(-1.0, -1.0);
    s_update_lnMolalityActCoeff();
    std::copy(m_lnActCoeffMolal_Unscaled.begin(), m_lnActCoeffMolal_Unscaled.end(), acMolality);
    for (size_t k = 0; k < m_kk; k++) {
        acMolality[k] = exp(acMolality[k]);
    }
}

// ------ Partial Molar Properties of the Solution -----------------

void HMWSoln::getChemPotentials(doublereal* mu) const
{
    double xx;

    // First get the standard chemical potentials in molar form. This requires
    // updates of standard state as a function of T and P
    getStandardChemPotentials(mu);

    // Update the activity coefficients. This also updates the internal molality
    // array.
    s_update_lnMolalityActCoeff();
    double xmolSolvent = moleFraction(0);
    for (size_t k = 1; k < m_kk; k++) {
        xx = std::max(m_molalities[k], SmallNumber);
        mu[k] += RT() * (log(xx) + m_lnActCoeffMolal_Scaled[k]);
    }
    xx = std::max(xmolSolvent, SmallNumber);
    mu[0] += RT() * (log(xx) + m_lnActCoeffMolal_Scaled[0]);
}

void HMWSoln::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // Get the nondimensional standard state enthalpies
    getEnthalpy_RT(hbar);

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RT() * temperature() * m_dlnActCoeffMolaldT_Scaled[k];
    }
}

void HMWSoln::getPartialMolarEntropies(doublereal* sbar) const
{
    // Get the standard state entropies at the temperature and pressure of the
    // solution.
    getEntropy_R(sbar);

    // Dimensionalize the entropies
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnMolalityActCoeff();

    // First we will add in the obvious dependence on the T term out front of
    // the log activity term
    doublereal mm;
    for (size_t k = 1; k < m_kk; k++) {
        mm = std::max(SmallNumber, m_molalities[k]);
        sbar[k] -= GasConstant * (log(mm) + m_lnActCoeffMolal_Scaled[k]);
    }
    double xmolSolvent = moleFraction(0);
    mm = std::max(SmallNumber, xmolSolvent);
    sbar[0] -= GasConstant *(log(mm) + m_lnActCoeffMolal_Scaled[0]);

    // Check to see whether activity coefficients are temperature dependent. If
    // they are, then calculate the their temperature derivatives and add them
    // into the result.
    s_update_dlnMolalityActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= RT() * m_dlnActCoeffMolaldT_Scaled[k];
    }
}

void HMWSoln::getPartialMolarVolumes(doublereal* vbar) const
{
    // Get the standard state values in m^3 kmol-1
    getStandardVolumes(vbar);

    // Update the derivatives wrt the activity coefficients.
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dP();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] += RT() * m_dlnActCoeffMolaldP_Scaled[k];
    }
}

void HMWSoln::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dT();
    s_update_d2lnMolalityActCoeff_dT2();
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] -= (2.0 * RT() * m_dlnActCoeffMolaldT_Scaled[k] +
                     RT() * temperature() * m_d2lnActCoeffMolaldT2_Scaled[k]);
    }
}

// -------------- Utilities -------------------------------

doublereal HMWSoln::satPressure(doublereal t) {
    double p_old = pressure();
    double t_old = temperature();
    double pres = m_waterSS->satPressure(t);

    // Set the underlying object back to its original state.
    m_waterSS->setState_TP(t_old, p_old);
    return pres;
}

static void check_nParams(const std::string& method, size_t nParams,
                          size_t m_formPitzerTemp)
{
    if (m_formPitzerTemp == PITZER_TEMP_CONSTANT && nParams != 1) {
        throw CanteraError(method, "'constant' temperature model requires one"
            " coefficient for each of parameter, but {} were given", nParams);
    } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR && nParams != 2) {
        throw CanteraError(method, "'linear' temperature model requires two"
            " coefficients for each parameter, but {} were given", nParams);
    }
    if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1 && nParams != 5) {
        throw CanteraError(method, "'complex' temperature model requires five"
            " coefficients for each parameter, but {} were given", nParams);
    }
}

void HMWSoln::setBinarySalt(const std::string& sp1, const std::string& sp2,
    size_t nParams, double* beta0, double* beta1, double* beta2,
    double* Cphi, double alpha1, double alpha2)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setBinarySalt", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setBinarySalt", "Species '{}' not found", sp2);
    }
    if (charge(k1) < 0 && charge(k2) > 0) {
        std::swap(k1, k2);
    } else if (charge(k1) * charge(k2) >= 0) {
        throw CanteraError("HMWSoln::setBinarySalt", "Species '{}' and '{}' "
            "do not have opposite charges ({}, {})", sp1, sp2,
            charge(k1), charge(k2));
    }
    check_nParams("HMWSoln::setBinarySalt", nParams, m_formPitzerTemp);

    size_t c = m_CounterIJ[k1 * m_kk + k2];
    m_Beta0MX_ij[c] = beta0[0];
    m_Beta1MX_ij[c] = beta1[0];
    m_Beta2MX_ij[c] = beta2[0];
    m_CphiMX_ij[c] = Cphi[0];
    for (size_t n = 0; n < nParams; n++) {
        m_Beta0MX_ij_coeff(n, c) = beta0[n];
        m_Beta1MX_ij_coeff(n, c) = beta1[n];
        m_Beta2MX_ij_coeff(n, c) = beta2[n];
        m_CphiMX_ij_coeff(n, c) = Cphi[n];
    }
    m_Alpha1MX_ij[c] = alpha1;
    m_Alpha2MX_ij[c] = alpha2;
}

void HMWSoln::setTheta(const std::string& sp1, const std::string& sp2,
        size_t nParams, double* theta)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setTheta", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setTheta", "Species '{}' not found", sp2);
    }
    if (charge(k1) * charge(k2) <= 0) {
        throw CanteraError("HMWSoln::setTheta", "Species '{}' and '{}' "
            "should both have the same (non-zero) charge ({}, {})", sp1, sp2,
            charge(k1), charge(k2));
    }
    check_nParams("HMWSoln::setTheta", nParams, m_formPitzerTemp);
    size_t c = m_CounterIJ[k1 * m_kk + k2];
    m_Theta_ij[c] = theta[0];
    for (size_t n = 0; n < nParams; n++) {
        m_Theta_ij_coeff(n, c) = theta[n];
    }
}

void HMWSoln::setPsi(const std::string& sp1, const std::string& sp2,
        const std::string& sp3, size_t nParams, double* psi)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    size_t k3 = speciesIndex(sp3);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setPsi", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setPsi", "Species '{}' not found", sp2);
    } else if (k3 == npos) {
        throw CanteraError("HMWSoln::setPsi", "Species '{}' not found", sp3);
    }

    if (!charge(k1) || !charge(k2) || !charge(k3) ||
        std::abs(sign(charge(k1) + sign(charge(k2)) + sign(charge(k3)))) != 1) {
        throw CanteraError("HMWSoln::setPsi", "All species must be ions and"
            " must include at least one cation and one anion, but given species"
            " (charges) were: {} ({}), {} ({}), and {} ({}).",
            sp1, charge(k1), sp2, charge(k2), sp3, charge(k3));
    }
    check_nParams("HMWSoln::setPsi", nParams, m_formPitzerTemp);
    auto cc = {k1*m_kk*m_kk + k2*m_kk + k3,
               k1*m_kk*m_kk + k3*m_kk + k2,
               k2*m_kk*m_kk + k1*m_kk + k3,
               k2*m_kk*m_kk + k3*m_kk + k1,
               k3*m_kk*m_kk + k2*m_kk + k1,
               k3*m_kk*m_kk + k1*m_kk + k2};
    for (auto c : cc) {
        for (size_t n = 0; n < nParams; n++) {
            m_Psi_ijk_coeff(n, c) = psi[n];
        }
        m_Psi_ijk[c] = psi[0];
    }
}

void HMWSoln::setLambda(const std::string& sp1, const std::string& sp2,
        size_t nParams, double* lambda)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setLambda", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setLambda", "Species '{}' not found", sp2);
    }

    if (charge(k1) != 0 && charge(k2) != 0) {
        throw CanteraError("HMWSoln::setLambda", "Expected at least one neutral"
            " species, but given species (charges) were: {} ({}) and {} ({}).",
            sp1, charge(k1), sp2, charge(k2));
    }
    if (charge(k1) != 0) {
        std::swap(k1, k2);
    }
    check_nParams("HMWSoln::setLambda", nParams, m_formPitzerTemp);
    size_t c = k1*m_kk + k2;
    for (size_t n = 0; n < nParams; n++) {
        m_Lambda_nj_coeff(n, c) = lambda[n];
    }
    m_Lambda_nj(k1, k2) = lambda[0];
}

void HMWSoln::setMunnn(const std::string& sp, size_t nParams, double* munnn)
{
    size_t k = speciesIndex(sp);
    if (k == npos) {
        throw CanteraError("HMWSoln::setMunnn", "Species '{}' not found", sp);
    }

    if (charge(k) != 0) {
        throw CanteraError("HMWSoln::setMunnn", "Expected a neutral species,"
                " got {} ({}).", sp, charge(k));
    }
    check_nParams("HMWSoln::setMunnn", nParams, m_formPitzerTemp);
    for (size_t n = 0; n < nParams; n++) {
        m_Mu_nnn_coeff(n, k) = munnn[n];
    }
    m_Mu_nnn[k] = munnn[0];
}

void HMWSoln::setZeta(const std::string& sp1, const std::string& sp2,
        const std::string& sp3, size_t nParams, double* psi)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    size_t k3 = speciesIndex(sp3);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setZeta", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setZeta", "Species '{}' not found", sp2);
    } else if (k3 == npos) {
        throw CanteraError("HMWSoln::setZeta", "Species '{}' not found", sp3);
    }

    if (charge(k1)*charge(k2)*charge(k3) != 0 ||
        sign(charge(k1)) + sign(charge(k2)) + sign(charge(k3)) != 0) {
        throw CanteraError("HMWSoln::setZeta", "Requires one neutral species, "
            "one cation, and one anion, but given species (charges) were: "
            "{} ({}), {} ({}), and {} ({}).",
            sp1, charge(k1), sp2, charge(k2), sp3, charge(k3));
    }

    //! Make k1 the neutral species
    if (charge(k2) == 0) {
        std::swap(k1, k2);
    } else if (charge(k3) == 0) {
        std::swap(k1, k3);
    }

    // Make k2 the cation
    if (charge(k3) > 0) {
        std::swap(k2, k3);
    }

    check_nParams("HMWSoln::setZeta", nParams, m_formPitzerTemp);
    // In contrast to setPsi, there are no duplicate entries
    size_t c = k1 * m_kk *m_kk + k2 * m_kk + k3;
    for (size_t n = 0; n < nParams; n++) {
        m_Psi_ijk_coeff(n, c) = psi[n];
    }
    m_Psi_ijk[c] = psi[0];
}

void HMWSoln::setPitzerTempModel(const std::string& model)
{
    if (caseInsensitiveEquals(model, "constant") || caseInsensitiveEquals(model, "default")) {
        m_formPitzerTemp = PITZER_TEMP_CONSTANT;
    } else if (caseInsensitiveEquals(model, "linear")) {
        m_formPitzerTemp = PITZER_TEMP_LINEAR;
    } else if (caseInsensitiveEquals(model, "complex") || caseInsensitiveEquals(model, "complex1")) {
        m_formPitzerTemp = PITZER_TEMP_COMPLEX1;
    } else {
        throw CanteraError("HMWSoln::setPitzerTempModel",
                           "Unknown Pitzer ActivityCoeff Temp model: {}", model);
    }
}

void HMWSoln::setA_Debye(double A)
{
    if (A < 0) {
        m_form_A_Debye = A_DEBYE_WATER;
    } else {
        m_form_A_Debye = A_DEBYE_CONST;
        m_A_Debye = A;
    }
}

void HMWSoln::setCroppingCoefficients(double ln_gamma_k_min,
    double ln_gamma_k_max, double ln_gamma_o_min, double ln_gamma_o_max)
{
        CROP_ln_gamma_k_min = ln_gamma_k_min;
        CROP_ln_gamma_k_max = ln_gamma_k_max;
        CROP_ln_gamma_o_min = ln_gamma_o_min;
        CROP_ln_gamma_o_max = ln_gamma_o_max;
}

vector_fp getSizedVector(const AnyMap& item, const std::string& key, size_t nCoeffs)
{
    vector_fp v;
    if (item[key].is<double>()) {
        // Allow a single value to be given directly, rather than as a list of
        // one item
        v.push_back(item[key].asDouble());
    } else {
        v = item[key].asVector<double>(1, nCoeffs);
    }
    if (v.size() == 1 && nCoeffs == 5) {
        // Adapt constant-temperature data to be compatible with the "complex"
        // temperature model
        v.resize(5, 0.0);
    }
    return v;
}

void HMWSoln::initThermo()
{
    MolalityVPSSTP::initThermo();
    if (m_input.hasKey("activity-data")) {
        auto& actData = m_input["activity-data"].as<AnyMap>();
        setPitzerTempModel(actData["temperature-model"].asString());
        initLengths();
        size_t nCoeffs = 1;
        if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
            nCoeffs = 2;
        } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
            nCoeffs = 5;
        }
        if (actData.hasKey("A_Debye")) {
            if (actData["A_Debye"] == "variable") {
                setA_Debye(-1);
            } else {
                setA_Debye(actData.convert("A_Debye", "kg^0.5/gmol^0.5"));
            }
        }
        if (actData.hasKey("max-ionic-strength")) {
            setMaxIonicStrength(actData["max-ionic-strength"].asDouble());
        }
        if (actData.hasKey("interactions")) {
            for (auto& item : actData["interactions"].asVector<AnyMap>()) {
                auto& species = item["species"].asVector<string>(1, 3);
                size_t nsp = species.size();
                double q0 = charge(speciesIndex(species[0]));
                double q1 = (nsp > 1) ? charge(speciesIndex(species[1])) : 0;
                double q2 = (nsp == 3) ? charge(speciesIndex(species[2])) : 0;
                if (nsp == 2 && q0 * q1 < 0) {
                    // Two species with opposite charges - binary salt
                    vector_fp beta0 = getSizedVector(item, "beta0", nCoeffs);
                    vector_fp beta1 = getSizedVector(item, "beta1", nCoeffs);
                    vector_fp beta2 = getSizedVector(item, "beta2", nCoeffs);
                    vector_fp Cphi = getSizedVector(item, "Cphi", nCoeffs);
                    if (beta0.size() != beta1.size() || beta0.size() != beta2.size()
                        || beta0.size() != Cphi.size()) {
                        throw InputFileError("HMWSoln::initThermo", item,
                            "Inconsistent binary salt array sizes ({}, {}, {}, {})",
                            beta0.size(), beta1.size(), beta2.size(), Cphi.size());
                    }
                    double alpha1 = item["alpha1"].asDouble();
                    double alpha2 = item.getDouble("alpha2", 0.0);
                    setBinarySalt(species[0], species[1], beta0.size(),
                        beta0.data(), beta1.data(), beta2.data(), Cphi.data(),
                        alpha1, alpha2);
                } else if (nsp == 2 && q0 * q1 > 0) {
                    // Two species with like charges - "theta" interaction
                    vector_fp theta = getSizedVector(item, "theta", nCoeffs);
                    setTheta(species[0], species[1], theta.size(), theta.data());
                } else if (nsp == 2 && q0 * q1 == 0) {
                    // Two species, including at least one neutral
                    vector_fp lambda = getSizedVector(item, "lambda", nCoeffs);
                    setLambda(species[0], species[1], lambda.size(), lambda.data());
                } else if (nsp == 3 && q0 * q1 * q2 != 0) {
                    // Three charged species - "psi" interaction
                    vector_fp psi = getSizedVector(item, "psi", nCoeffs);
                    setPsi(species[0], species[1], species[2],
                           psi.size(), psi.data());
                } else if (nsp == 3 && q0 * q1 * q2 == 0) {
                    // Three species, including one neutral
                    vector_fp zeta = getSizedVector(item, "zeta", nCoeffs);
                    setZeta(species[0], species[1], species[2],
                            zeta.size(), zeta.data());
                } else if (nsp == 1) {
                    // single species (should be neutral)
                    vector_fp mu = getSizedVector(item, "mu", nCoeffs);
                    setMunnn(species[0], mu.size(), mu.data());
                }
            }
        }
        if (actData.hasKey("cropping-coefficients")) {
            auto& crop = actData["cropping-coefficients"].as<AnyMap>();
            setCroppingCoefficients(
                crop.getDouble("ln_gamma_k_min", -5.0),
                crop.getDouble("ln_gamma_k_max", 15.0),
                crop.getDouble("ln_gamma_o_min", -6.0),
                crop.getDouble("ln_gamma_o_max", 3.0));
        }
    } else {
        initLengths();
    }

    for (int i = 0; i < 17; i++) {
        elambda[i] = 0.0;
        elambda1[i] = 0.0;
    }

    // Store a local pointer to the water standard state model.
    m_waterSS = providePDSS(0);

    // Initialize the water property calculator. It will share the internal eos
    // water calculator.
    m_waterProps.reset(new WaterProps(dynamic_cast<PDSS_Water*>(m_waterSS)));

    // Lastly calculate the charge balance and then add stuff until the charges
    // compensate
    vector_fp mf(m_kk, 0.0);
    getMoleFractions(mf.data());
    bool notDone = true;

    while (notDone) {
        double sum = 0.0;
        size_t kMaxC = npos;
        double MaxC = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            sum += mf[k] * charge(k);
            if (fabs(mf[k] * charge(k)) > MaxC) {
                kMaxC = k;
            }
        }
        size_t kHp = speciesIndex("H+");
        size_t kOHm = speciesIndex("OH-");

        if (fabs(sum) > 1.0E-30) {
            if (kHp != npos) {
                if (mf[kHp] > sum * 1.1) {
                    mf[kHp] -= sum;
                    mf[0] += sum;
                    notDone = false;
                } else {
                    if (sum > 0.0) {
                        mf[kHp] *= 0.5;
                        mf[0] += mf[kHp];
                        sum -= mf[kHp];
                    }
                }
            }
            if (notDone) {
                if (kOHm != npos) {
                    if (mf[kOHm] > -sum * 1.1) {
                        mf[kOHm] += sum;
                        mf[0] -= sum;
                        notDone = false;
                    } else {
                        if (sum < 0.0) {
                            mf[kOHm] *= 0.5;
                            mf[0] += mf[kOHm];
                            sum += mf[kOHm];
                        }
                    }
                }
                if (notDone && kMaxC != npos) {
                    if (mf[kMaxC] > (1.1 * sum / charge(kMaxC))) {
                        mf[kMaxC] -= sum / charge(kMaxC);
                        mf[0] += sum / charge(kMaxC);
                    } else {
                        mf[kMaxC] *= 0.5;
                        mf[0] += mf[kMaxC];
                        notDone = true;
                    }
                }
            }
            setMoleFractions(mf.data());
        } else {
            notDone = false;
        }
    }

    calcIMSCutoffParams_();
    calcMCCutoffParams_();
    setMoleFSolventMin(1.0E-5);
}

void HMWSoln::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id_) {
            throw CanteraError("HMWSoln::initThermoXML",
                               "phasenode and Id are incompatible");
        }
    }

    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("HMWSoln::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    // Determine the form of the Pitzer model, We will use this information to
    // size arrays below.
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& scNode = thermoNode.child("activityCoefficients");

        // Determine the form of the temperature dependence of the Pitzer
        // activity coefficient model.
        string formString = scNode.attrib("TempModel");
        if (formString != "") {
            setPitzerTempModel(formString);
        }

        // Determine the reference temperature of the Pitzer activity
        // coefficient model's temperature dependence formulation: defaults to
        // 25C
        formString = scNode.attrib("TempReference");
        if (formString != "") {
            setPitzerRefTemperature(fpValueCheck(formString));
        }
    }

    // Initialize all of the lengths of arrays in the object
    // now that we know what species are in the phase.
    initLengths();

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");

        // Look for parameters for A_Debye
        if (acNode.hasChild("A_Debye")) {
            XML_Node& ADebye = acNode.child("A_Debye");
            if (caseInsensitiveEquals(ADebye["model"], "water")) {
                setA_Debye(-1);
            } else {
                setA_Debye(getFloat(acNode, "A_Debye"));
            }
        }

        // Look for Parameters for the Maximum Ionic Strength
        if (acNode.hasChild("maxIonicStrength")) {
            setMaxIonicStrength(getFloat(acNode, "maxIonicStrength"));
        }

        for (const auto& xmlACChild : acNode.children()) {
            string nodeName = xmlACChild->name();

            // Process any of the XML fields that make up the Pitzer Database.
            // Entries will be ignored if any of the species in the entry aren't
            // in the solution.
            if (caseInsensitiveEquals(nodeName, "binarysaltparameters")) {
                readXMLBinarySalt(*xmlACChild);
            } else if (caseInsensitiveEquals(nodeName, "thetaanion")) {
                readXMLTheta(*xmlACChild);
            } else if (caseInsensitiveEquals(nodeName, "thetacation")) {
                readXMLTheta(*xmlACChild);
            } else if (caseInsensitiveEquals(nodeName, "psicommonanion")) {
                readXMLPsi(*xmlACChild);
            } else if (caseInsensitiveEquals(nodeName, "psicommoncation")) {
                readXMLPsi(*xmlACChild);
            } else if (caseInsensitiveEquals(nodeName, "lambdaneutral")) {
                readXMLLambdaNeutral(*xmlACChild);
            } else if (caseInsensitiveEquals(nodeName, "zetacation")) {
                readXMLZetaCation(*xmlACChild);
            }
        }

        // Go look up the optional Cropping parameters
        if (acNode.hasChild("croppingCoefficients")) {
            XML_Node& cropNode = acNode.child("croppingCoefficients");
            setCroppingCoefficients(
                getFloat(cropNode.child("ln_gamma_k_min"), "pureSolventValue"),
                getFloat(cropNode.child("ln_gamma_k_max"), "pureSolventValue"),
                getFloat(cropNode.child("ln_gamma_o_min"), "pureSolventValue"),
                getFloat(cropNode.child("ln_gamma_o_max"), "pureSolventValue"));
        }
    }

    MolalityVPSSTP::initThermoXML(phaseNode, id_);
}

double HMWSoln::A_Debye_TP(double tempArg, double presArg) const
{
    double T = temperature();
    double A;
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }

    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if(cached.validate(T, P)) {
        return m_A_Debye;
    }

    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        A = m_A_Debye;
        break;
    case A_DEBYE_WATER:
        A = m_waterProps->ADebye(T, P, 0);
        m_A_Debye = A;
        break;
    default:
        throw CanteraError("HMWSoln::A_Debye_TP", "shouldn't be here");
    }
    return A;
}

double HMWSoln::dA_DebyedT_TP(double tempArg, double presArg) const
{
    doublereal T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    doublereal P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    doublereal dAdT;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdT = 0.0;
        break;
    case A_DEBYE_WATER:
        dAdT = m_waterProps->ADebye(T, P, 1);
        break;
    default:
        throw CanteraError("HMWSoln::dA_DebyedT_TP", "shouldn't be here");
    }
    return dAdT;
}

double HMWSoln::dA_DebyedP_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }

    double dAdP;
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdP = 0.0;
        break;
    case A_DEBYE_WATER:
        if(cached.validate(T, P)) {
            dAdP = cached.value;
        } else {
            dAdP = m_waterProps->ADebye(T, P, 3);
            cached.value = dAdP;
        }
        break;
    default:
        throw CanteraError("HMWSoln::dA_DebyedP_TP", "shouldn't be here");
    }
    return dAdP;
}

double HMWSoln::ADebye_L(double tempArg, double presArg) const
{
    double dAdT = dA_DebyedT_TP();
    double dAphidT = dAdT /3.0;
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    return dAphidT * (4.0 * GasConstant * T * T);
}

double HMWSoln::ADebye_V(double tempArg, double presArg) const
{
    double dAdP = dA_DebyedP_TP();
    double dAphidP = dAdP /3.0;
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    return - dAphidP * (4.0 * GasConstant * T);
}

double HMWSoln::ADebye_J(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double A_L = ADebye_L(T, presArg);
    double d2 = d2A_DebyedT2_TP(T, presArg);
    double d2Aphi = d2 / 3.0;
    return 2.0 * A_L / T + 4.0 * GasConstant * T * T *d2Aphi;
}

double HMWSoln::d2A_DebyedT2_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double d2AdT2;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        d2AdT2 = 0.0;
        break;
    case A_DEBYE_WATER:
        d2AdT2 = m_waterProps->ADebye(T, P, 2);
        break;
    default:
        throw CanteraError("HMWSoln::d2A_DebyedT2_TP", "shouldn't be here");
    }
    return d2AdT2;
}

// ---------- Other Property Functions

// ------------ Private and Restricted Functions ------------------

void HMWSoln::initLengths()
{
    m_tmpV.resize(m_kk, 0.0);
    m_molalitiesCropped.resize(m_kk, 0.0);

    size_t maxCounterIJlen = 1 + (m_kk-1) * (m_kk-2) / 2;

    // Figure out the size of the temperature coefficient arrays
    int TCoeffLength = 1;
    if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
        TCoeffLength = 2;
    } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        TCoeffLength = 5;
    }

    m_Beta0MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Beta1MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Beta2MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_CphiMX_ij.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_L.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_P.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Alpha1MX_ij.resize(maxCounterIJlen, 2.0);
    m_Alpha2MX_ij.resize(maxCounterIJlen, 12.0);
    m_Theta_ij.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_L.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_P.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    size_t n = m_kk*m_kk*m_kk;
    m_Psi_ijk.resize(n, 0.0);
    m_Psi_ijk_L.resize(n, 0.0);
    m_Psi_ijk_LL.resize(n, 0.0);
    m_Psi_ijk_P.resize(n, 0.0);
    m_Psi_ijk_coeff.resize(TCoeffLength, n, 0.0);

    m_Lambda_nj.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_L.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_LL.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_P.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_coeff.resize(TCoeffLength, m_kk * m_kk, 0.0);

    m_Mu_nnn.resize(m_kk, 0.0);
    m_Mu_nnn_L.resize(m_kk, 0.0);
    m_Mu_nnn_LL.resize(m_kk, 0.0);
    m_Mu_nnn_P.resize(m_kk, 0.0);
    m_Mu_nnn_coeff.resize(TCoeffLength, m_kk, 0.0);

    m_lnActCoeffMolal_Scaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldT_Scaled.resize(m_kk, 0.0);
    m_d2lnActCoeffMolaldT2_Scaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldP_Scaled.resize(m_kk, 0.0);
    m_lnActCoeffMolal_Unscaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldT_Unscaled.resize(m_kk, 0.0);
    m_d2lnActCoeffMolaldT2_Unscaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldP_Unscaled.resize(m_kk, 0.0);

    m_CounterIJ.resize(m_kk*m_kk, 0);
    m_gfunc_IJ.resize(maxCounterIJlen, 0.0);
    m_g2func_IJ.resize(maxCounterIJlen, 0.0);
    m_hfunc_IJ.resize(maxCounterIJlen, 0.0);
    m_h2func_IJ.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_L.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_P.resize(maxCounterIJlen, 0.0);
    m_Phiprime_IJ.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_L.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_P.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_P.resize(maxCounterIJlen, 0.0);

    m_gamma_tmp.resize(m_kk, 0.0);
    IMS_lnActCoeffMolal_.resize(m_kk, 0.0);
    CROP_speciesCropped_.resize(m_kk, 0);

    counterIJ_setup();
}

void HMWSoln::s_update_lnMolalityActCoeff() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    // Calculate the molalities. Currently, the molalities may not be current
    // with respect to the contents of the State objects' data.
    calcMolalities();

    // Calculate a cropped set of molalities that will be used in all activity
    // coefficient calculations.
    calcMolalitiesCropped();

    // Update the temperature dependence of the Pitzer coefficients and their
    // derivatives
    s_updatePitzer_CoeffWRTemp();

    // Calculate the IMS cutoff factors
    s_updateIMS_lnMolalityActCoeff();

    // Now do the main calculation.
    s_updatePitzer_lnMolalityActCoeff();

    double xmolSolvent = moleFraction(0);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);
    double lnActCoeffMolal0 = - log(xx) + (xx - 1.0)/xx;
    double lnxs = log(xx);

    for (size_t k = 1; k < m_kk; k++) {
        CROP_speciesCropped_[k] = 0;
        m_lnActCoeffMolal_Unscaled[k] += IMS_lnActCoeffMolal_[k];
        if (m_lnActCoeffMolal_Unscaled[k] > (CROP_ln_gamma_k_max- 2.5 *lnxs)) {
            CROP_speciesCropped_[k] = 2;
            m_lnActCoeffMolal_Unscaled[k] = CROP_ln_gamma_k_max - 2.5 * lnxs;
        }
        if (m_lnActCoeffMolal_Unscaled[k] < (CROP_ln_gamma_k_min - 2.5 *lnxs)) {
            // -1.0 and -1.5 caused multiple solutions
            CROP_speciesCropped_[k] = 2;
            m_lnActCoeffMolal_Unscaled[k] = CROP_ln_gamma_k_min - 2.5 * lnxs;
        }
    }
    CROP_speciesCropped_[0] = 0;
    m_lnActCoeffMolal_Unscaled[0] += (IMS_lnActCoeffMolal_[0] - lnActCoeffMolal0);
    if (m_lnActCoeffMolal_Unscaled[0] < CROP_ln_gamma_o_min) {
        CROP_speciesCropped_[0] = 2;
        m_lnActCoeffMolal_Unscaled[0] = CROP_ln_gamma_o_min;
    }
    if (m_lnActCoeffMolal_Unscaled[0] > CROP_ln_gamma_o_max) {
        CROP_speciesCropped_[0] = 2;
        // -0.5 caused multiple solutions
        m_lnActCoeffMolal_Unscaled[0] = CROP_ln_gamma_o_max;
    }
    if (m_lnActCoeffMolal_Unscaled[0] > CROP_ln_gamma_o_max - 0.5 * lnxs) {
        CROP_speciesCropped_[0] = 2;
        m_lnActCoeffMolal_Unscaled[0] = CROP_ln_gamma_o_max - 0.5 * lnxs;
    }

    // Now do the pH Scaling
    s_updateScaling_pHScaling();
}

void HMWSoln::calcMolalitiesCropped() const
{
    doublereal Imax = 0.0;
    m_molalitiesAreCropped = false;

    for (size_t k = 0; k < m_kk; k++) {
        m_molalitiesCropped[k] = m_molalities[k];
        Imax = std::max(m_molalities[k] * charge(k) * charge(k), Imax);
    }

    int cropMethod = 1;
    if (cropMethod == 0) {
        // Quick return
        if (Imax < m_maxIionicStrength) {
            return;
        }

        m_molalitiesAreCropped = true;
        for (size_t i = 1; i < (m_kk - 1); i++) {
            double charge_i = charge(i);
            double abs_charge_i = fabs(charge_i);
            if (charge_i == 0.0) {
                continue;
            }
            for (size_t j = (i+1); j < m_kk; j++) {
                double charge_j = charge(j);
                double abs_charge_j = fabs(charge_j);

                // Only loop over oppositely charge species
                if (charge_i * charge_j < 0) {
                    double Iac_max = m_maxIionicStrength;

                    if (m_molalitiesCropped[i] > m_molalitiesCropped[j]) {
                        Imax = m_molalitiesCropped[i] * abs_charge_i * abs_charge_i;
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[i] = Iac_max / (abs_charge_i * abs_charge_i);
                        }
                        Imax = m_molalitiesCropped[j] * fabs(abs_charge_j * abs_charge_i);
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[j] = Iac_max / (abs_charge_j * abs_charge_i);
                        }
                    } else {
                        Imax = m_molalitiesCropped[j] * abs_charge_j * abs_charge_j;
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[j] = Iac_max / (abs_charge_j * abs_charge_j);
                        }
                        Imax = m_molalitiesCropped[i] * abs_charge_j * abs_charge_i;
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[i] = Iac_max / (abs_charge_j * abs_charge_i);
                        }
                    }
                }
            }
        }

        // Do this loop 10 times until we have achieved charge neutrality in the
        // cropped molalities
        for (int times = 0; times< 10; times++) {
            double anion_charge = 0.0;
            double cation_charge = 0.0;
            size_t anion_contrib_max_i = npos;
            double anion_contrib_max = -1.0;
            size_t cation_contrib_max_i = npos;
            double cation_contrib_max = -1.0;
            for (size_t i = 0; i < m_kk; i++) {
                double charge_i = charge(i);
                if (charge_i < 0.0) {
                    double anion_contrib = - m_molalitiesCropped[i] * charge_i;
                    anion_charge += anion_contrib;
                    if (anion_contrib > anion_contrib_max) {
                        anion_contrib_max = anion_contrib;
                        anion_contrib_max_i = i;
                    }
                } else if (charge_i > 0.0) {
                    double cation_contrib = m_molalitiesCropped[i] * charge_i;
                    cation_charge += cation_contrib;
                    if (cation_contrib > cation_contrib_max) {
                        cation_contrib_max = cation_contrib;
                        cation_contrib_max_i = i;
                    }
                }
            }
            double total_charge = cation_charge - anion_charge;
            if (total_charge > 1.0E-8) {
                double desiredCrop = total_charge/charge(cation_contrib_max_i);
                double maxCrop = 0.66 * m_molalitiesCropped[cation_contrib_max_i];
                if (desiredCrop < maxCrop) {
                    m_molalitiesCropped[cation_contrib_max_i] -= desiredCrop;
                    break;
                } else {
                    m_molalitiesCropped[cation_contrib_max_i] -= maxCrop;
                }
            } else if (total_charge < -1.0E-8) {
                double desiredCrop = total_charge/charge(anion_contrib_max_i);
                double maxCrop = 0.66 * m_molalitiesCropped[anion_contrib_max_i];
                if (desiredCrop < maxCrop) {
                    m_molalitiesCropped[anion_contrib_max_i] -= desiredCrop;
                    break;
                } else {
                    m_molalitiesCropped[anion_contrib_max_i] -= maxCrop;
                }
            } else {
                break;
            }
        }
    }

    if (cropMethod == 1) {
        double* molF = m_gamma_tmp.data();
        getMoleFractions(molF);
        double xmolSolvent = molF[0];
        if (xmolSolvent >= MC_X_o_cutoff_) {
            return;
        }

        m_molalitiesAreCropped = true;
        double poly = MC_apCut_ + MC_bpCut_ * xmolSolvent + MC_dpCut_* xmolSolvent * xmolSolvent;
        double p = xmolSolvent + MC_epCut_ + exp(- xmolSolvent/ MC_cpCut_) * poly;
        double denomInv = 1.0/ (m_Mnaught * p);
        for (size_t k = 0; k < m_kk; k++) {
            m_molalitiesCropped[k] = molF[k] * denomInv;
        }

        // Do a further check to see if the Ionic strength is below a max value
        // Reduce the molalities to enforce this. Note, this algorithm preserves
        // the charge neutrality of the solution after cropping.
        double Itmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            Itmp += m_molalitiesCropped[k] * charge(k) * charge(k);
        }
        if (Itmp > m_maxIionicStrength) {
            double ratio = Itmp / m_maxIionicStrength;
            for (size_t k = 0; k < m_kk; k++) {
                if (charge(k) != 0.0) {
                    m_molalitiesCropped[k] *= ratio;
                }
            }
        }
    }
}

void HMWSoln::counterIJ_setup() const
{
    m_CounterIJ.resize(m_kk * m_kk);
    int counter = 0;
    for (size_t i = 0; i < m_kk; i++) {
        size_t n = i;
        size_t nc = m_kk * i;
        m_CounterIJ[n] = 0;
        m_CounterIJ[nc] = 0;
    }
    for (size_t i = 1; i < (m_kk - 1); i++) {
        size_t n = m_kk * i + i;
        m_CounterIJ[n] = 0;
        for (size_t j = (i+1); j < m_kk; j++) {
            n = m_kk * j + i;
            size_t nc = m_kk * i + j;
            counter++;
            m_CounterIJ[n] = counter;
            m_CounterIJ[nc] = counter;
        }
    }
}

void HMWSoln::readXMLBinarySalt(XML_Node& BinSalt)
{
    if (BinSalt.name() != "binarySaltParameters") {
        throw CanteraError("HMWSoln::readXMLBinarySalt",
                           "Incorrect name for processing this routine: " + BinSalt.name());
    }

    string iName = BinSalt.attrib("cation");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "no cation attrib");
    }
    string jName = BinSalt.attrib("anion");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "no anion attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos) {
        return;
    }

    vector_fp beta0, beta1, beta2, Cphi;
    getFloatArray(BinSalt, beta0, false, "", "beta0");
    getFloatArray(BinSalt, beta1, false, "", "beta1");
    getFloatArray(BinSalt, beta2, false, "", "beta2");
    getFloatArray(BinSalt, Cphi, false, "", "Cphi");
    if (beta0.size() != beta1.size() || beta0.size() != beta2.size() ||
        beta0.size() != Cphi.size()) {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "Inconsistent"
            " array sizes ({}, {}, {}, {})", beta0.size(), beta1.size(),
            beta2.size(), Cphi.size());
    }
    if (beta0.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        beta0.resize(5, 0.0);
        beta1.resize(5, 0.0);
        beta2.resize(5, 0.0);
        Cphi.resize(5, 0.0);
    }
    double alpha1 = getFloat(BinSalt, "Alpha1");
    double alpha2 = 0.0;
    getOptionalFloat(BinSalt, "Alpha2", alpha2);
    setBinarySalt(iName, jName, beta0.size(), beta0.data(), beta1.data(),
        beta2.data(), Cphi.data(), alpha1, alpha2);
}

void HMWSoln::readXMLTheta(XML_Node& node)
{
    string ispName, jspName;
    if (node.name() == "thetaAnion") {
        ispName = node.attrib("anion1");
        if (ispName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no anion1 attrib");
        }
        jspName = node.attrib("anion2");
        if (jspName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no anion2 attrib");
        }
    } else if (node.name() == "thetaCation") {
        ispName = node.attrib("cation1");
        if (ispName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no cation1 attrib");
        }
        jspName = node.attrib("cation2");
        if (jspName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no cation2 attrib");
        }
    } else {
        throw CanteraError("HMWSoln::readXMLTheta",
                "Incorrect name for processing this routine: " + node.name());
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(ispName) == npos || speciesIndex(jspName) == npos) {
        return;
    }

    vector_fp theta;
    getFloatArray(node, theta, false, "", "theta");
    if (theta.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        theta.resize(5, 0.0);
    }
    setTheta(ispName, jspName, theta.size(), theta.data());
}

void HMWSoln::readXMLPsi(XML_Node& node)
{
    string iName, jName, kName;
    if (node.name() == "psiCommonCation") {
        kName = node.attrib("cation");
        if (kName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no cation attrib");
        }
        iName = node.attrib("anion1");
        if (iName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no anion1 attrib");
        }
        jName = node.attrib("anion2");
        if (jName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no anion2 attrib");
        }
    } else if (node.name() == "psiCommonAnion") {
        kName = node.attrib("anion");
        if (kName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no anion attrib");
        }
        iName = node.attrib("cation1");
        if (iName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no cation1 attrib");
        }
        jName = node.attrib("cation2");
        if (jName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no cation2 attrib");
        }
    } else {
        throw CanteraError("HMWSoln::readXMLPsi",
                   "Incorrect name for processing this routine: " + node.name());
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos ||
        speciesIndex(kName) == npos) {
        return;
    }

    vector_fp psi;
    getFloatArray(node, psi, false, "", "psi");
    if (psi.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        psi.resize(5, 0.0);
    }
    setPsi(iName, jName, kName, psi.size(), psi.data());
}

void HMWSoln::readXMLLambdaNeutral(XML_Node& node)
{
    vector_fp vParams;
    if (node.name() != "lambdaNeutral") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral",
                           "Incorrect name for processing this routine: " + node.name());
    }
    string iName = node.attrib("species1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no species1 attrib");
    }
    string jName = node.attrib("species2");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no species2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos) {
        return;
    }

    vector_fp lambda;
    getFloatArray(node, lambda, false, "", "lambda");
    if (lambda.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        lambda.resize(5, 0.0);
    }
    setLambda(iName, jName, lambda.size(), lambda.data());
}

void HMWSoln::readXMLMunnnNeutral(XML_Node& node)
{
    if (node.name() != "MunnnNeutral") {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral",
                           "Incorrect name for processing this routine: " + node.name());
    }
    string iName = node.attrib("species1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral", "no species1 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos) {
        return;
    }

    vector_fp munnn;
    getFloatArray(node, munnn, false, "", "munnn");
    if (munnn.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        munnn.resize(5, 0.0);
    }
    setMunnn(iName, munnn.size(), munnn.data());
}

void HMWSoln::readXMLZetaCation(const XML_Node& node)
{
    if (node.name() != "zetaCation") {
        throw CanteraError("HMWSoln::readXMLZetaCation",
                           "Incorrect name for processing this routine: " + node.name());
    }

    string iName = node.attrib("neutral");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no neutral attrib");
    }
    string jName = node.attrib("cation1");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no cation1 attrib");
    }
    string kName = node.attrib("anion1");
    if (kName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no anion1 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos ||
        speciesIndex(kName) == npos) {
        return;
    }

    vector_fp zeta;
    getFloatArray(node, zeta, false, "", "zeta");
    if (zeta.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        zeta.resize(5, 0.0);
    }
    setZeta(iName, jName, kName, zeta.size(), zeta.data());
}

void HMWSoln::calcIMSCutoffParams_()
{
    double IMS_gamma_o_min_ = 1.0E-5; // value at the zero solvent point
    double IMS_gamma_k_min_ = 10.0; // minimum at the zero solvent point
    double IMS_slopefCut_ = 0.6; // slope of the f function at the zero solvent point

    IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_);
    IMS_efCut_ = 0.0;
    bool converged = false;
    double oldV = 0.0;
    for (int its = 0; its < 100 && !converged; its++) {
        oldV = IMS_efCut_;
        IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_) -IMS_efCut_;
        IMS_bfCut_ = IMS_afCut_ / IMS_cCut_ + IMS_slopefCut_ - 1.0;
        IMS_dfCut_ = ((- IMS_afCut_/IMS_cCut_ + IMS_bfCut_ - IMS_bfCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_afCut_ + IMS_X_o_cutoff_*(IMS_bfCut_ + IMS_dfCut_ *IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_efCut_ = - eterm * tmp;
        if (fabs(IMS_efCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("HMWSoln::calcIMSCutoffParams_()",
                           " failed to converge on the f polynomial");
    }
    converged = false;
    double f_0 = IMS_afCut_ + IMS_efCut_;
    double f_prime_0 = 1.0 - IMS_afCut_ / IMS_cCut_ + IMS_bfCut_;
    IMS_egCut_ = 0.0;
    for (int its = 0; its < 100 && !converged; its++) {
        oldV = IMS_egCut_;
        double lng_0 = -log(IMS_gamma_o_min_) - f_prime_0 / f_0;
        IMS_agCut_ = exp(lng_0) - IMS_egCut_;
        IMS_bgCut_ = IMS_agCut_ / IMS_cCut_ + IMS_slopegCut_ - 1.0;
        IMS_dgCut_ = ((- IMS_agCut_/IMS_cCut_ + IMS_bgCut_ - IMS_bgCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_agCut_ + IMS_X_o_cutoff_*(IMS_bgCut_ + IMS_dgCut_ *IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_egCut_ = - eterm * tmp;
        if (fabs(IMS_egCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("HMWSoln::calcIMSCutoffParams_()",
                           " failed to converge on the g polynomial");
    }
}

void HMWSoln::calcMCCutoffParams_()
{
    double MC_X_o_min_ = 0.35; // value at the zero solvent point
    MC_X_o_cutoff_ = 0.6;
    double MC_slopepCut_ = 0.02; // slope of the p function at the zero solvent point
    MC_cpCut_ = 0.25;

    // Initial starting values
    MC_apCut_ = MC_X_o_min_;
    MC_epCut_ = 0.0;
    bool converged = false;
    double oldV = 0.0;
    double damp = 0.5;
    for (int its = 0; its < 500 && !converged; its++) {
        oldV = MC_epCut_;
        MC_apCut_ = damp *(MC_X_o_min_ - MC_epCut_) + (1-damp) * MC_apCut_;
        double MC_bpCutNew = MC_apCut_ / MC_cpCut_ + MC_slopepCut_ - 1.0;
        MC_bpCut_ = damp * MC_bpCutNew + (1-damp) * MC_bpCut_;
        double MC_dpCutNew = ((- MC_apCut_/MC_cpCut_ + MC_bpCut_ - MC_bpCut_ * MC_X_o_cutoff_/MC_cpCut_)
                              /
                              (MC_X_o_cutoff_ * MC_X_o_cutoff_/MC_cpCut_ - 2.0 * MC_X_o_cutoff_));
        MC_dpCut_ = damp * MC_dpCutNew + (1-damp) * MC_dpCut_;
        double tmp = MC_apCut_ + MC_X_o_cutoff_*(MC_bpCut_ + MC_dpCut_ * MC_X_o_cutoff_);
        double eterm = std::exp(- MC_X_o_cutoff_ / MC_cpCut_);
        MC_epCut_ = - eterm * tmp;
        double diff = MC_epCut_ - oldV;
        if (fabs(diff) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("HMWSoln::calcMCCutoffParams_()",
                           " failed to converge on the p polynomial");
    }
}

void HMWSoln::s_updatePitzer_CoeffWRTemp(int doDerivs) const
{
    double T = temperature();
    const double twoT = 2.0 * T;
    const double invT = 1.0 / T;
    const double invT2 = invT * invT;
    const double twoinvT3 = 2.0 * invT * invT2;
    double tinv = 0.0, tln = 0.0, tlin = 0.0, tquad = 0.0;
    if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
        tlin = T - m_TempPitzerRef;
    } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        tlin = T - m_TempPitzerRef;
        tquad = T * T - m_TempPitzerRef * m_TempPitzerRef;
        tln = log(T/ m_TempPitzerRef);
        tinv = 1.0/T - 1.0/m_TempPitzerRef;
    }

    for (size_t i = 1; i < (m_kk - 1); i++) {
        for (size_t j = (i+1); j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            const double* beta0MX_coeff = m_Beta0MX_ij_coeff.ptrColumn(counterIJ);
            const double* beta1MX_coeff = m_Beta1MX_ij_coeff.ptrColumn(counterIJ);
            const double* beta2MX_coeff = m_Beta2MX_ij_coeff.ptrColumn(counterIJ);
            const double* CphiMX_coeff = m_CphiMX_ij_coeff.ptrColumn(counterIJ);
            const double* Theta_coeff = m_Theta_ij_coeff.ptrColumn(counterIJ);

            switch (m_formPitzerTemp) {
            case PITZER_TEMP_CONSTANT:
                break;
            case PITZER_TEMP_LINEAR:

                m_Beta0MX_ij[counterIJ] = beta0MX_coeff[0]
                                          + beta0MX_coeff[1]*tlin;
                m_Beta0MX_ij_L[counterIJ] = beta0MX_coeff[1];
                m_Beta0MX_ij_LL[counterIJ] = 0.0;
                m_Beta1MX_ij[counterIJ] = beta1MX_coeff[0]
                                            + beta1MX_coeff[1]*tlin;
                m_Beta1MX_ij_L[counterIJ] = beta1MX_coeff[1];
                m_Beta1MX_ij_LL[counterIJ] = 0.0;
                m_Beta2MX_ij[counterIJ] = beta2MX_coeff[0]
                                             + beta2MX_coeff[1]*tlin;
                m_Beta2MX_ij_L[counterIJ] = beta2MX_coeff[1];
                m_Beta2MX_ij_LL[counterIJ] = 0.0;
                m_CphiMX_ij[counterIJ] = CphiMX_coeff[0]
                                             + CphiMX_coeff[1]*tlin;
                m_CphiMX_ij_L[counterIJ] = CphiMX_coeff[1];
                m_CphiMX_ij_LL[counterIJ] = 0.0;
                m_Theta_ij[counterIJ] = Theta_coeff[0] + Theta_coeff[1]*tlin;
                m_Theta_ij_L[counterIJ] = Theta_coeff[1];
                m_Theta_ij_LL[counterIJ] = 0.0;
                break;

            case PITZER_TEMP_COMPLEX1:
                m_Beta0MX_ij[counterIJ] = beta0MX_coeff[0]
                                          + beta0MX_coeff[1]*tlin
                                          + beta0MX_coeff[2]*tquad
                                          + beta0MX_coeff[3]*tinv
                                          + beta0MX_coeff[4]*tln;
                m_Beta1MX_ij[counterIJ] = beta1MX_coeff[0]
                                          + beta1MX_coeff[1]*tlin
                                          + beta1MX_coeff[2]*tquad
                                          + beta1MX_coeff[3]*tinv
                                          + beta1MX_coeff[4]*tln;
                m_Beta2MX_ij[counterIJ] = beta2MX_coeff[0]
                                          + beta2MX_coeff[1]*tlin
                                          + beta2MX_coeff[2]*tquad
                                          + beta2MX_coeff[3]*tinv
                                          + beta2MX_coeff[4]*tln;
                m_CphiMX_ij[counterIJ] = CphiMX_coeff[0]
                                         + CphiMX_coeff[1]*tlin
                                         + CphiMX_coeff[2]*tquad
                                         + CphiMX_coeff[3]*tinv
                                         + CphiMX_coeff[4]*tln;
                m_Theta_ij[counterIJ] = Theta_coeff[0]
                                        + Theta_coeff[1]*tlin
                                        + Theta_coeff[2]*tquad
                                        + Theta_coeff[3]*tinv
                                        + Theta_coeff[4]*tln;
                m_Beta0MX_ij_L[counterIJ] = beta0MX_coeff[1]
                                             + beta0MX_coeff[2]*twoT
                                             - beta0MX_coeff[3]*invT2
                                             + beta0MX_coeff[4]*invT;
                m_Beta1MX_ij_L[counterIJ] = beta1MX_coeff[1]
                                             + beta1MX_coeff[2]*twoT
                                             - beta1MX_coeff[3]*invT2
                                             + beta1MX_coeff[4]*invT;
                m_Beta2MX_ij_L[counterIJ] = beta2MX_coeff[1]
                                             + beta2MX_coeff[2]*twoT
                                             - beta2MX_coeff[3]*invT2
                                             + beta2MX_coeff[4]*invT;
                m_CphiMX_ij_L[counterIJ] = CphiMX_coeff[1]
                                            + CphiMX_coeff[2]*twoT
                                            - CphiMX_coeff[3]*invT2
                                            + CphiMX_coeff[4]*invT;
                m_Theta_ij_L[counterIJ] = Theta_coeff[1]
                                            + Theta_coeff[2]*twoT
                                            - Theta_coeff[3]*invT2
                                            + Theta_coeff[4]*invT;
                doDerivs = 2;
                if (doDerivs > 1) {
                    m_Beta0MX_ij_LL[counterIJ] =
                        + beta0MX_coeff[2]*2.0
                        + beta0MX_coeff[3]*twoinvT3
                        - beta0MX_coeff[4]*invT2;
                    m_Beta1MX_ij_LL[counterIJ] =
                        + beta1MX_coeff[2]*2.0
                        + beta1MX_coeff[3]*twoinvT3
                        - beta1MX_coeff[4]*invT2;
                    m_Beta2MX_ij_LL[counterIJ] =
                        + beta2MX_coeff[2]*2.0
                        + beta2MX_coeff[3]*twoinvT3
                        - beta2MX_coeff[4]*invT2;
                    m_CphiMX_ij_LL[counterIJ] =
                        + CphiMX_coeff[2]*2.0
                        + CphiMX_coeff[3]*twoinvT3
                        - CphiMX_coeff[4]*invT2;
                    m_Theta_ij_LL[counterIJ] =
                        + Theta_coeff[2]*2.0
                        + Theta_coeff[3]*twoinvT3
                        - Theta_coeff[4]*invT2;
                }
                break;
            }
        }
    }

    // Lambda interactions and Mu_nnn
    // i must be neutral for this term to be nonzero. We take advantage of this
    // here to lower the operation count.
    for (size_t i = 1; i < m_kk; i++) {
        if (charge(i) == 0.0) {
            for (size_t j = 1; j < m_kk; j++) {
                size_t n = i * m_kk + j;
                const double* Lambda_coeff = m_Lambda_nj_coeff.ptrColumn(n);
                switch (m_formPitzerTemp) {
                case PITZER_TEMP_CONSTANT:
                    m_Lambda_nj(i,j) = Lambda_coeff[0];
                    break;
                case PITZER_TEMP_LINEAR:
                    m_Lambda_nj(i,j) = Lambda_coeff[0] + Lambda_coeff[1]*tlin;
                    m_Lambda_nj_L(i,j) = Lambda_coeff[1];
                    m_Lambda_nj_LL(i,j) = 0.0;
                    break;
                case PITZER_TEMP_COMPLEX1:
                    m_Lambda_nj(i,j) = Lambda_coeff[0]
                                       + Lambda_coeff[1]*tlin
                                       + Lambda_coeff[2]*tquad
                                       + Lambda_coeff[3]*tinv
                                       + Lambda_coeff[4]*tln;

                    m_Lambda_nj_L(i,j) = Lambda_coeff[1]
                                         + Lambda_coeff[2]*twoT
                                         - Lambda_coeff[3]*invT2
                                         + Lambda_coeff[4]*invT;

                    m_Lambda_nj_LL(i,j) =
                        Lambda_coeff[2]*2.0
                        + Lambda_coeff[3]*twoinvT3
                        - Lambda_coeff[4]*invT2;
                }

                if (j == i) {
                    const double* Mu_coeff = m_Mu_nnn_coeff.ptrColumn(i);
                    switch (m_formPitzerTemp) {
                    case PITZER_TEMP_CONSTANT:
                        m_Mu_nnn[i] = Mu_coeff[0];
                        break;
                    case PITZER_TEMP_LINEAR:
                        m_Mu_nnn[i] = Mu_coeff[0] + Mu_coeff[1]*tlin;
                        m_Mu_nnn_L[i] = Mu_coeff[1];
                        m_Mu_nnn_LL[i] = 0.0;
                        break;
                    case PITZER_TEMP_COMPLEX1:
                        m_Mu_nnn[i] = Mu_coeff[0]
                                      + Mu_coeff[1]*tlin
                                      + Mu_coeff[2]*tquad
                                      + Mu_coeff[3]*tinv
                                      + Mu_coeff[4]*tln;
                        m_Mu_nnn_L[i] = Mu_coeff[1]
                                        + Mu_coeff[2]*twoT
                                        - Mu_coeff[3]*invT2
                                        + Mu_coeff[4]*invT;
                        m_Mu_nnn_LL[i] =
                            Mu_coeff[2]*2.0
                            + Mu_coeff[3]*twoinvT3
                            - Mu_coeff[4]*invT2;
                    }
                }
            }
        }
    }

    switch(m_formPitzerTemp) {
    case PITZER_TEMP_CONSTANT:
      for (size_t i = 1; i < m_kk; i++) {
          for (size_t j = 1; j < m_kk; j++) {
              for (size_t k = 1; k < m_kk; k++) {
                  size_t n = i * m_kk *m_kk + j * m_kk + k;
                  const double* Psi_coeff = m_Psi_ijk_coeff.ptrColumn(n);
                  m_Psi_ijk[n] = Psi_coeff[0];
              }
          }
      }
      break;
    case PITZER_TEMP_LINEAR:
      for (size_t i = 1; i < m_kk; i++) {
          for (size_t j = 1; j < m_kk; j++) {
              for (size_t k = 1; k < m_kk; k++) {
                  size_t n = i * m_kk *m_kk + j * m_kk + k;
                  const double* Psi_coeff = m_Psi_ijk_coeff.ptrColumn(n);
                  m_Psi_ijk[n] = Psi_coeff[0] + Psi_coeff[1]*tlin;
                  m_Psi_ijk_L[n] = Psi_coeff[1];
                  m_Psi_ijk_LL[n] = 0.0;
              }
          }
      }
      break;
    case PITZER_TEMP_COMPLEX1:
      for (size_t i = 1; i < m_kk; i++) {
          for (size_t j = 1; j < m_kk; j++) {
              for (size_t k = 1; k < m_kk; k++) {
                  size_t n = i * m_kk *m_kk + j * m_kk + k;
                  const double* Psi_coeff = m_Psi_ijk_coeff.ptrColumn(n);
                  m_Psi_ijk[n] = Psi_coeff[0]
                                 + Psi_coeff[1]*tlin
                                 + Psi_coeff[2]*tquad
                                 + Psi_coeff[3]*tinv
                                 + Psi_coeff[4]*tln;
                  m_Psi_ijk_L[n] = Psi_coeff[1]
                                   + Psi_coeff[2]*twoT
                                   - Psi_coeff[3]*invT2
                                   + Psi_coeff[4]*invT;
                  m_Psi_ijk_LL[n] =
                      Psi_coeff[2]*2.0
                      + Psi_coeff[3]*twoinvT3
                      - Psi_coeff[4]*invT2;
              }
          }
      }
      break;
    }
}

void HMWSoln::s_updatePitzer_lnMolalityActCoeff() const
{
    // Use the CROPPED molality of the species in solution.
    const vector_fp& molality = m_molalitiesCropped;

    // These are data inputs about the Pitzer correlation. They come from the
    // input file for the Pitzer model.
    vector_fp& gamma_Unscaled = m_gamma_tmp;

    // Local variables defined by Coltrin
    double etheta[5][5], etheta_prime[5][5], sqrtIs;

    // Molality based ionic strength of the solution
    double Is = 0.0;

    // Molar charge of the solution: In Pitzer's notation, this is his variable
    // called "Z".
    double molarcharge = 0.0;

    // molalitysum is the sum of the molalities over all solutes, even those
    // with zero charge.
    double molalitysumUncropped = 0.0;

    // Make sure the counter variables are setup
    counterIJ_setup();

    // ---------- Calculate common sums over solutes ---------------------
    for (size_t n = 1; n < m_kk; n++) {
        // ionic strength
        Is += charge(n) * charge(n) * molality[n];
        // total molar charge
        molarcharge += fabs(charge(n)) * molality[n];
        molalitysumUncropped += m_molalities[n];
    }
    Is *= 0.5;

    // Store the ionic molality in the object for reference.
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);

    // The following call to calc_lambdas() calculates all 16 elements of the
    // elambda and elambda1 arrays, given the value of the ionic strength (Is)
    calc_lambdas(Is);

    // Step 2: Find the coefficients E-theta and E-thetaprime for all
    // combinations of positive unlike charges up to 4
    for (int z1 = 1; z1 <=4; z1++) {
        for (int z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
        }
    }

    // calculate g(x) and hfunc(x) for each cation-anion pair MX. In the
    // original literature, hfunc, was called gprime. However, it's not the
    // derivative of g(x), so I renamed it.
    for (size_t i = 1; i < (m_kk - 1); i++) {
        for (size_t j = (i+1); j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // Only loop over oppositely charge species
            if (charge(i)*charge(j) < 0) {
                // x is a reduced function variable
                double x1 = sqrtIs * m_Alpha1MX_ij[counterIJ];
                if (x1 > 1.0E-100) {
                    m_gfunc_IJ[counterIJ] = 2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 * x1);
                    m_hfunc_IJ[counterIJ] = -2.0 *
                                       (1.0-(1.0 + x1 + 0.5 * x1 * x1) * exp(-x1)) / (x1 * x1);
                } else {
                    m_gfunc_IJ[counterIJ] = 0.0;
                    m_hfunc_IJ[counterIJ] = 0.0;
                }

                if (m_Beta2MX_ij[counterIJ] != 0.0) {
                    double x2 = sqrtIs * m_Alpha2MX_ij[counterIJ];
                    if (x2 > 1.0E-100) {
                        m_g2func_IJ[counterIJ] = 2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        m_h2func_IJ[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        m_g2func_IJ[counterIJ] = 0.0;
                        m_h2func_IJ[counterIJ] = 0.0;
                    }
                }
            } else {
                m_gfunc_IJ[counterIJ] = 0.0;
                m_hfunc_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION TO CALCULATE BMX, BprimeMX, BphiMX
    // Agrees with Pitzer, Eq. (49), (51), (55)
    for (size_t i = 1; i < m_kk - 1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive and the
            // other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_BMX_IJ[counterIJ] = m_Beta0MX_ij[counterIJ]
                                  + m_Beta1MX_ij[counterIJ] * m_gfunc_IJ[counterIJ]
                                  + m_Beta2MX_ij[counterIJ] * m_g2func_IJ[counterIJ];

                if (Is > 1.0E-150) {
                    m_BprimeMX_IJ[counterIJ] = (m_Beta1MX_ij[counterIJ] * m_hfunc_IJ[counterIJ]/Is +
                                           m_Beta2MX_ij[counterIJ] * m_h2func_IJ[counterIJ]/Is);
                } else {
                    m_BprimeMX_IJ[counterIJ] = 0.0;
                }
                m_BphiMX_IJ[counterIJ] = m_BMX_IJ[counterIJ] + Is*m_BprimeMX_IJ[counterIJ];
            } else {
                m_BMX_IJ[counterIJ] = 0.0;
                m_BprimeMX_IJ[counterIJ] = 0.0;
                m_BphiMX_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION TO CALCULATE CMX
    // Agrees with Pitzer, Eq. (53).
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_CMX_IJ[counterIJ] = m_CphiMX_ij[counterIJ]/
                                 (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                m_CMX_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi
    // Agrees with Pitzer, Eq. 72, 73, 74
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive and the
            // other is negative
            if (charge(i)*charge(j) > 0) {
                int z1 = (int) fabs(charge(i));
                int z2 = (int) fabs(charge(j));
                m_Phi_IJ[counterIJ] = m_Theta_ij[counterIJ] + etheta[z1][z2];
                m_Phiprime_IJ[counterIJ] = etheta_prime[z1][z2];
                m_PhiPhi_IJ[counterIJ] = m_Phi_IJ[counterIJ] + Is * m_Phiprime_IJ[counterIJ];
            } else {
                m_Phi_IJ[counterIJ] = 0.0;
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION FOR CALCULATION OF F
    // Agrees with Pitzer Eqn. (65)
    double Aphi = A_Debye_TP() / 3.0;
    double F = -Aphi * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                 + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive and the
            // other is negative
            if (charge(i)*charge(j) < 0) {
                F += molality[i]*molality[j] * m_BprimeMX_IJ[counterIJ];
            }

            // Both species have a non-zero charge, and they
            // have the same sign
            if (charge(i)*charge(j) > 0) {
                F += molality[i]*molality[j] * m_Phiprime_IJ[counterIJ];
            }
        }
    }

    for (size_t i = 1; i < m_kk; i++) {

        // SUBSECTION FOR CALCULATING THE ACTCOEFF FOR CATIONS
        // equations agree with my notes, Eqn. (118).
        // Equations agree with Pitzer, eqn.(63)
        if (charge(i) > 0.0) {
            // species i is the cation (positive) to calc the actcoeff
            double zsqF = charge(i)*charge(i)*F;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 += molality[j] *
                            (2.0*m_BMX_IJ[counterIJ] + molarcharge*m_CMX_IJ[counterIJ]);
                    if (j < m_kk-1) {
                        // This term is the ternary interaction involving the
                        // non-duplicate sum over double anions, j, k, with
                        // respect to the cation, i.
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk[n];
                            }
                        }
                    }
                }

                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk[n];

                            // Find the counterIJ for the j,k interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += (fabs(charge(i))*
                                           molality[j]*molality[k]*m_CMX_IJ[counterIJ2]);
                        }
                    }
                }

                // Handle neutral j species
                if (charge(j) == 0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj(j,i);

                    // Zeta interaction term
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t izeta = j;
                            size_t jzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                            double zeta = m_Psi_ijk[n];
                            if (zeta != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta;
                            }
                        }
                    }
                }
            }

            // Add all of the contributions up to yield the log of the solute
            // activity coefficients (molality scale)
            m_lnActCoeffMolal_Unscaled[i] = zsqF + sum1 + sum2 + sum3 + sum4 + sum5;
            gamma_Unscaled[i] = exp(m_lnActCoeffMolal_Unscaled[i]);
        }

        // SUBSECTION FOR CALCULATING THE ACTCOEFF FOR ANIONS
        // equations agree with my notes, Eqn. (119).
        // Equations agree with Pitzer, eqn.(64)
        if (charge(i) < 0) {
            // species i is an anion (negative)
            double zsqF = charge(i)*charge(i)*F;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                // For Anions, do the cation interactions.
                if (charge(j) > 0) {
                    sum1 += molality[j]*
                           (2.0*m_BMX_IJ[counterIJ]+molarcharge*m_CMX_IJ[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk[n];
                            }
                        }
                    }
                }

                // For Anions, do the other anion interactions.
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk[n];
                            // Find the counterIJ for the symmetric binary interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i))*
                                    molality[j]*molality[k]*m_CMX_IJ[counterIJ2];
                        }
                    }
                }

                // for Anions, do the neutral species interaction
                if (charge(j) == 0.0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj(j,i);
                    // Zeta interaction term
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta = m_Psi_ijk[n];
                            if (zeta != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta;
                            }
                        }
                    }
                }
            }
            m_lnActCoeffMolal_Unscaled[i] = zsqF + sum1 + sum2 + sum3 + sum4 + sum5;
            gamma_Unscaled[i] = exp(m_lnActCoeffMolal_Unscaled[i]);
        }

        // SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF
        // equations agree with my notes,
        // Equations agree with Pitzer,
        if (charge(i) == 0.0) {
            double sum1 = 0.0;
            double sum3 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                sum1 += molality[j]*2.0*m_Lambda_nj(i,j);
                // Zeta term -> we piggyback on the psi term
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*m_Psi_ijk[n];
                        }
                    }
                }
            }
            double sum2 = 3.0 * molality[i]* molality[i] * m_Mu_nnn[i];
            m_lnActCoeffMolal_Unscaled[i] = sum1 + sum2 + sum3;
            gamma_Unscaled[i] = exp(m_lnActCoeffMolal_Unscaled[i]);
        }
    }

    // SUBSECTION FOR CALCULATING THE OSMOTIC COEFF
    // equations agree with my notes, Eqn. (117).
    // Equations agree with Pitzer, eqn.(62)
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;

    // term1 is the DH term in the osmotic coefficient expression
    // b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer
    //                          implementations.
    // Is = Ionic strength on the molality scale (units of (gmol/kg))
    // Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
    double term1 = -Aphi * pow(Is,1.5) / (1.0 + 1.2 * sqrt(Is));

    for (size_t j = 1; j < m_kk; j++) {
        // Loop Over Cations
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];

                    sum1 += molality[j]*molality[k]*
                            (m_BphiMX_IJ[counterIJ] + molarcharge*m_CMX_IJ[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_lnMolalityActCoeff",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between 2 cations.
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum2 += molality[j]*molality[k]*m_PhiPhi_IJ[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*molality[m]*m_Psi_ijk[n];
                        }
                    }
                }
            }
        }

        // Loop Over Anions
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_lnMolalityActCoeff",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between two anions
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum3 += molality[j]*molality[k]*m_PhiPhi_IJ[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*molality[m]*m_Psi_ijk[n];
                        }
                    }
                }
            }
        }

        // Loop Over Neutral Species
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 += molality[j]*molality[k]*m_Lambda_nj(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 += molality[j]*molality[k]*m_Lambda_nj(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 += molality[j]*molality[k]*m_Lambda_nj(j,k);
                    } else if (k == j) {
                        sum6 += 0.5 * molality[j]*molality[k]*m_Lambda_nj(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            size_t n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta = m_Psi_ijk[n];
                            if (zeta != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta;
                            }
                        }
                    }
                }
            }
            sum7 += molality[j]*molality[j]*molality[j]*m_Mu_nnn[j];
        }
    }
    double sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);
    // Calculate the osmotic coefficient from
    //     osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
    double osmotic_coef;
    if (molalitysumUncropped > 1.0E-150) {
        osmotic_coef = 1.0 + (sum_m_phi_minus_1 / molalitysumUncropped);
    } else {
        osmotic_coef = 1.0;
    }
    double lnwateract = -(m_weightSolvent/1000.0) * molalitysumUncropped * osmotic_coef;

    // In Cantera, we define the activity coefficient of the solvent as
    //
    //     act_0 = actcoeff_0 * Xmol_0
    //
    // We have just computed act_0. However, this routine returns
    //     ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
    double xmolSolvent = moleFraction(0);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);
    m_lnActCoeffMolal_Unscaled[0] = lnwateract - log(xx);
}

void HMWSoln::s_update_dlnMolalityActCoeff_dT() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    // Zero the unscaled 2nd derivatives
    m_dlnActCoeffMolaldT_Unscaled.assign(m_kk, 0.0);

    // Do the actual calculation of the unscaled temperature derivatives
    s_updatePitzer_dlnMolalityActCoeff_dT();

    for (size_t k = 1; k < m_kk; k++) {
        if (CROP_speciesCropped_[k] == 2) {
            m_dlnActCoeffMolaldT_Unscaled[k] = 0.0;
        }
    }

    if (CROP_speciesCropped_[0]) {
        m_dlnActCoeffMolaldT_Unscaled[0] = 0.0;
    }

    // Do the pH scaling to the derivatives
    s_updateScaling_pHScaling_dT();
}

void HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT() const
{
    // It may be assumed that the Pitzer activity coefficient routine is called
    // immediately preceding the calling of this routine. Therefore, some
    // quantities do not need to be recalculated in this routine.

    const vector_fp& molality = m_molalitiesCropped;
    double* d_gamma_dT_Unscaled = m_gamma_tmp.data();

    // Local variables defined by Coltrin
    double etheta[5][5], etheta_prime[5][5], sqrtIs;

    // Molality based ionic strength of the solution
    double Is = 0.0;

    // Molar charge of the solution: In Pitzer's notation, this is his variable
    // called "Z".
    double molarcharge = 0.0;

    // molalitysum is the sum of the molalities over all solutes, even those
    // with zero charge.
    double molalitysum = 0.0;

    // Make sure the counter variables are setup
    counterIJ_setup();

    // ---------- Calculate common sums over solutes ---------------------
    for (size_t n = 1; n < m_kk; n++) {
        // ionic strength
        Is += charge(n) * charge(n) * molality[n];
        // total molar charge
        molarcharge += fabs(charge(n)) * molality[n];
        molalitysum += molality[n];
    }
    Is *= 0.5;

    // Store the ionic molality in the object for reference.
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);

    // The following call to calc_lambdas() calculates all 16 elements of the
    // elambda and elambda1 arrays, given the value of the ionic strength (Is)
    calc_lambdas(Is);

    // Step 2: Find the coefficients E-theta and E-thetaprime for all
    // combinations of positive unlike charges up to 4
    for (int z1 = 1; z1 <=4; z1++) {
        for (int z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
        }
    }

    // calculate g(x) and hfunc(x) for each cation-anion pair MX
    // In the original literature, hfunc, was called gprime. However,
    // it's not the derivative of g(x), so I renamed it.
    for (size_t i = 1; i < (m_kk - 1); i++) {
        for (size_t j = (i+1); j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // Only loop over oppositely charge species
            if (charge(i)*charge(j) < 0) {
                // x is a reduced function variable
                double x1 = sqrtIs * m_Alpha1MX_ij[counterIJ];
                if (x1 > 1.0E-100) {
                    m_gfunc_IJ[counterIJ] = 2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 * x1);
                    m_hfunc_IJ[counterIJ] = -2.0 *
                                       (1.0-(1.0 + x1 + 0.5 * x1 *x1) * exp(-x1)) / (x1 * x1);
                } else {
                    m_gfunc_IJ[counterIJ] = 0.0;
                    m_hfunc_IJ[counterIJ] = 0.0;
                }

                if (m_Beta2MX_ij_L[counterIJ] != 0.0) {
                    double x2 = sqrtIs * m_Alpha2MX_ij[counterIJ];
                    if (x2 > 1.0E-100) {
                        m_g2func_IJ[counterIJ] = 2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        m_h2func_IJ[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        m_g2func_IJ[counterIJ] = 0.0;
                        m_h2func_IJ[counterIJ] = 0.0;
                    }
                }
            } else {
                m_gfunc_IJ[counterIJ] = 0.0;
                m_hfunc_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION TO CALCULATE BMX_L, BprimeMX_L, BphiMX_L
    // These are now temperature derivatives of the previously calculated
    // quantities.
    for (size_t i = 1; i < m_kk - 1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_BMX_IJ_L[counterIJ] = m_Beta0MX_ij_L[counterIJ]
                                    + m_Beta1MX_ij_L[counterIJ] * m_gfunc_IJ[counterIJ]
                                    + m_Beta2MX_ij_L[counterIJ] * m_gfunc_IJ[counterIJ];
                if (Is > 1.0E-150) {
                    m_BprimeMX_IJ_L[counterIJ] = (m_Beta1MX_ij_L[counterIJ] * m_hfunc_IJ[counterIJ]/Is +
                                             m_Beta2MX_ij_L[counterIJ] * m_h2func_IJ[counterIJ]/Is);
                } else {
                    m_BprimeMX_IJ_L[counterIJ] = 0.0;
                }
                m_BphiMX_IJ_L[counterIJ] = m_BMX_IJ_L[counterIJ] + Is*m_BprimeMX_IJ_L[counterIJ];
            } else {
                m_BMX_IJ_L[counterIJ] = 0.0;
                m_BprimeMX_IJ_L[counterIJ] = 0.0;
                m_BphiMX_IJ_L[counterIJ] = 0.0;
            }
        }
    }

    // --------- SUBSECTION TO CALCULATE CMX_L ----------
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_CMX_IJ_L[counterIJ] = m_CphiMX_ij_L[counterIJ]/
                                   (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                m_CMX_IJ_L[counterIJ] = 0.0;
            }
        }
    }

    // ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) > 0) {
                m_Phi_IJ_L[counterIJ] = m_Theta_ij_L[counterIJ];
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ_L[counterIJ] = m_Phi_IJ_L[counterIJ] + Is * m_Phiprime_IJ[counterIJ];
            } else {
                m_Phi_IJ_L[counterIJ] = 0.0;
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ_L[counterIJ] = 0.0;
            }
        }
    }

    // ----------- SUBSECTION FOR CALCULATION OF dFdT ---------------------
    double dA_DebyedT = dA_DebyedT_TP();
    double dAphidT = dA_DebyedT /3.0;
    double dFdT = -dAphidT * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                       + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0) {
                dFdT += molality[i]*molality[j] * m_BprimeMX_IJ_L[counterIJ];
            }

            // Both species have a non-zero charge, and they
            // have the same sign, e.g., both positive or both negative.
            if (charge(i)*charge(j) > 0) {
                dFdT += molality[i]*molality[j] * m_Phiprime_IJ[counterIJ];
            }
        }
    }

    for (size_t i = 1; i < m_kk; i++) {
        // -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
        if (charge(i) > 0) {
            // species i is the cation (positive) to calc the actcoeff
            double zsqdFdT = charge(i)*charge(i)*dFdT;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 += molality[j]*
                            (2.0*m_BMX_IJ_L[counterIJ] + molarcharge*m_CMX_IJ_L[counterIJ]);
                    if (j < m_kk-1) {
                        // This term is the ternary interaction involving the
                        // non-duplicate sum over double anions, j, k, with
                        // respect to the cation, i.
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk_L[n];
                            }
                        }
                    }
                }

                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ_L[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk_L[n];

                            // Find the counterIJ for the j,k interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i))*
                                    molality[j]*molality[k]*m_CMX_IJ_L[counterIJ2];
                        }
                    }
                }

                // Handle neutral j species
                if (charge(j) == 0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj_L(j,i);
                }

                // Zeta interaction term
                for (size_t k = 1; k < m_kk; k++) {
                    if (charge(k) < 0.0) {
                        size_t izeta = j;
                        size_t jzeta = i;
                        n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                        double zeta_L = m_Psi_ijk_L[n];
                        if (zeta_L != 0.0) {
                            sum5 += molality[j]*molality[k]*zeta_L;
                        }
                    }
                }
            }

            // Add all of the contributions up to yield the log of the
            // solute activity coefficients (molality scale)
            m_dlnActCoeffMolaldT_Unscaled[i] =
                zsqdFdT + sum1 + sum2 + sum3 + sum4 + sum5;
            d_gamma_dT_Unscaled[i] = exp(m_dlnActCoeffMolaldT_Unscaled[i]);
        }

        // ------ SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR ANIONS ------
        if (charge(i) < 0) {
            // species i is an anion (negative)
            double zsqdFdT = charge(i)*charge(i)*dFdT;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                // For Anions, do the cation interactions.
                if (charge(j) > 0) {
                    sum1 += molality[j]*
                            (2.0*m_BMX_IJ_L[counterIJ] + molarcharge*m_CMX_IJ_L[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk_L[n];
                            }
                        }
                    }
                }

                // For Anions, do the other anion interactions.
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ_L[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk_L[n];
                            // Find the counterIJ for the symmetric binary interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i)) *
                                    molality[j]*molality[k]*m_CMX_IJ_L[counterIJ2];
                        }
                    }
                }

                // for Anions, do the neutral species interaction
                if (charge(j) == 0.0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj_L(j,i);
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta_L = m_Psi_ijk_L[n];
                            if (zeta_L != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta_L;
                            }
                        }
                    }
                }
            }
            m_dlnActCoeffMolaldT_Unscaled[i] =
                zsqdFdT + sum1 + sum2 + sum3 + sum4 + sum5;
            d_gamma_dT_Unscaled[i] = exp(m_dlnActCoeffMolaldT_Unscaled[i]);
        }

        // SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF
        // equations agree with my notes,
        // Equations agree with Pitzer,
        if (charge(i) == 0.0) {
            double sum1 = 0.0;
            double sum3 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                sum1 += molality[j]*2.0*m_Lambda_nj_L(i,j);
                // Zeta term -> we piggyback on the psi term
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*m_Psi_ijk_L[n];
                        }
                    }
                }
            }
            double sum2 = 3.0 * molality[i] * molality[i] * m_Mu_nnn_L[i];
            m_dlnActCoeffMolaldT_Unscaled[i] = sum1 + sum2 + sum3;
            d_gamma_dT_Unscaled[i] = exp(m_dlnActCoeffMolaldT_Unscaled[i]);
        }
    }

    // ------ SUBSECTION FOR CALCULATING THE d OSMOTIC COEFF dT ---------
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;

    // term1 is the temperature derivative of the DH term in the osmotic
    // coefficient expression
    // b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer implementations.
    // Is = Ionic strength on the molality scale (units of (gmol/kg))
    // Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
    double term1 = -dAphidT * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (size_t j = 1; j < m_kk; j++) {
        // Loop Over Cations
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum1 += molality[j]*molality[k]*
                            (m_BphiMX_IJ_L[counterIJ] + molarcharge*m_CMX_IJ_L[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between 2 cations.
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum2 += molality[j]*molality[k]*m_PhiPhi_IJ_L[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*molality[m]*m_Psi_ijk_L[n];
                        }
                    }
                }
            }
        }

        // Loop Over Anions
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between two anions
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum3 += molality[j]*molality[k]*m_PhiPhi_IJ_L[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*molality[m]*m_Psi_ijk_L[n];
                        }
                    }
                }
            }
        }

        // Loop Over Neutral Species
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 += molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 += molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 += molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                    } else if (k == j) {
                        sum6 += 0.5 * molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            size_t n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta_L = m_Psi_ijk_L[n];
                            if (zeta_L != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta_L;
                            }
                        }
                    }
                }
            }
            sum7 += molality[j]*molality[j]*molality[j]*m_Mu_nnn_L[j];
        }
    }
    double sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);
    // Calculate the osmotic coefficient from
    //     osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
    double d_osmotic_coef_dT;
    if (molalitysum > 1.0E-150) {
        d_osmotic_coef_dT = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
        d_osmotic_coef_dT = 0.0;
    }

    double d_lnwateract_dT = -(m_weightSolvent/1000.0) * molalitysum * d_osmotic_coef_dT;

    // In Cantera, we define the activity coefficient of the solvent as
    //
    //     act_0 = actcoeff_0 * Xmol_0
    //
    // We have just computed act_0. However, this routine returns
    //     ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
    m_dlnActCoeffMolaldT_Unscaled[0] = d_lnwateract_dT;
}

void HMWSoln::s_update_d2lnMolalityActCoeff_dT2() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    // Zero the unscaled 2nd derivatives
    m_d2lnActCoeffMolaldT2_Unscaled.assign(m_kk, 0.0);

    //! Calculate the unscaled 2nd derivatives
    s_updatePitzer_d2lnMolalityActCoeff_dT2();

    for (size_t k = 1; k < m_kk; k++) {
        if (CROP_speciesCropped_[k] == 2) {
            m_d2lnActCoeffMolaldT2_Unscaled[k] = 0.0;
        }
    }

    if (CROP_speciesCropped_[0]) {
        m_d2lnActCoeffMolaldT2_Unscaled[0] = 0.0;
    }

    // Scale the 2nd derivatives
    s_updateScaling_pHScaling_dT2();
}

void HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2() const
{
    const double* molality = m_molalitiesCropped.data();

    // Local variables defined by Coltrin
    double etheta[5][5], etheta_prime[5][5], sqrtIs;

    // Molality based ionic strength of the solution
    double Is = 0.0;

    // Molar charge of the solution: In Pitzer's notation, this is his variable
    // called "Z".
    double molarcharge = 0.0;

    // molalitysum is the sum of the molalities over all solutes, even those
    // with zero charge.
    double molalitysum = 0.0;

    // Make sure the counter variables are setup
    counterIJ_setup();

    // ---------- Calculate common sums over solutes ---------------------
    for (size_t n = 1; n < m_kk; n++) {
        // ionic strength
        Is += charge(n) * charge(n) * molality[n];
        // total molar charge
        molarcharge += fabs(charge(n)) * molality[n];
        molalitysum += molality[n];
    }
    Is *= 0.5;

    // Store the ionic molality in the object for reference.
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);

    // The following call to calc_lambdas() calculates all 16 elements of the
    // elambda and elambda1 arrays, given the value of the ionic strength (Is)
    calc_lambdas(Is);

    // Step 2: Find the coefficients E-theta and E-thetaprime for all
    // combinations of positive unlike charges up to 4
    for (int z1 = 1; z1 <=4; z1++) {
        for (int z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
        }
    }

    // calculate gfunc(x) and hfunc(x) for each cation-anion pair MX. In the
    // original literature, hfunc, was called gprime. However, it's not the
    // derivative of gfunc(x), so I renamed it.
    for (size_t i = 1; i < (m_kk - 1); i++) {
        for (size_t j = (i+1); j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // Only loop over oppositely charge species
            if (charge(i)*charge(j) < 0) {
                // x is a reduced function variable
                double x1 = sqrtIs * m_Alpha1MX_ij[counterIJ];
                if (x1 > 1.0E-100) {
                    m_gfunc_IJ[counterIJ] = 2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 *x1);
                    m_hfunc_IJ[counterIJ] = -2.0*
                                       (1.0-(1.0 + x1 + 0.5*x1 * x1) * exp(-x1)) / (x1 * x1);
                } else {
                    m_gfunc_IJ[counterIJ] = 0.0;
                    m_hfunc_IJ[counterIJ] = 0.0;
                }

                if (m_Beta2MX_ij_LL[counterIJ] != 0.0) {
                    double x2 = sqrtIs * m_Alpha2MX_ij[counterIJ];
                    if (x2 > 1.0E-100) {
                        m_g2func_IJ[counterIJ] = 2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        m_h2func_IJ[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        m_g2func_IJ[counterIJ] = 0.0;
                        m_h2func_IJ[counterIJ] = 0.0;
                    }
                }
            } else {
                m_gfunc_IJ[counterIJ] = 0.0;
                m_hfunc_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION TO CALCULATE BMX_L, BprimeMX_LL, BphiMX_L
    // These are now temperature derivatives of the previously calculated
    // quantities.
    for (size_t i = 1; i < m_kk - 1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_BMX_IJ_LL[counterIJ] = m_Beta0MX_ij_LL[counterIJ]
                                     + m_Beta1MX_ij_LL[counterIJ] * m_gfunc_IJ[counterIJ]
                                     + m_Beta2MX_ij_LL[counterIJ] * m_g2func_IJ[counterIJ];
                if (Is > 1.0E-150) {
                    m_BprimeMX_IJ_LL[counterIJ] = (m_Beta1MX_ij_LL[counterIJ] * m_hfunc_IJ[counterIJ]/Is +
                                              m_Beta2MX_ij_LL[counterIJ] * m_h2func_IJ[counterIJ]/Is);
                } else {
                    m_BprimeMX_IJ_LL[counterIJ] = 0.0;
                }
                m_BphiMX_IJ_LL[counterIJ] = m_BMX_IJ_LL[counterIJ] + Is*m_BprimeMX_IJ_LL[counterIJ];
            } else {
                m_BMX_IJ_LL[counterIJ] = 0.0;
                m_BprimeMX_IJ_LL[counterIJ] = 0.0;
                m_BphiMX_IJ_LL[counterIJ] = 0.0;
            }
        }
    }

    // --------- SUBSECTION TO CALCULATE CMX_LL ----------
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_CMX_IJ_LL[counterIJ] = m_CphiMX_ij_LL[counterIJ]/
                                    (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                m_CMX_IJ_LL[counterIJ] = 0.0;
            }
        }
    }

    // ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) > 0) {
                m_Phi_IJ_LL[counterIJ] = m_Theta_ij_LL[counterIJ];
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ_LL[counterIJ] = m_Phi_IJ_LL[counterIJ];
            } else {
                m_Phi_IJ_LL[counterIJ] = 0.0;
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ_LL[counterIJ] = 0.0;
            }
        }
    }

    // ----------- SUBSECTION FOR CALCULATION OF d2FdT2 ---------------------
    double d2AphidT2 = d2A_DebyedT2_TP() / 3.0;
    double d2FdT2 = -d2AphidT2 * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                           + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0) {
                d2FdT2 += molality[i]*molality[j] * m_BprimeMX_IJ_LL[counterIJ];
            }

            // Both species have a non-zero charge, and they
            // have the same sign, e.g., both positive or both negative.
            if (charge(i)*charge(j) > 0) {
                d2FdT2 += molality[i]*molality[j] * m_Phiprime_IJ[counterIJ];
            }
        }
    }

    for (size_t i = 1; i < m_kk; i++) {
        // -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
        if (charge(i) > 0) {
            // species i is the cation (positive) to calc the actcoeff
            double zsqd2FdT2 = charge(i)*charge(i)*d2FdT2;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 += molality[j]*
                            (2.0*m_BMX_IJ_LL[counterIJ] + molarcharge*m_CMX_IJ_LL[counterIJ]);
                    if (j < m_kk-1) {
                        // This term is the ternary interaction involving the
                        // non-duplicate sum over double anions, j, k, with
                        // respect to the cation, i.
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk_LL[n];
                            }
                        }
                    }
                }

                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ_LL[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk_LL[n];

                            // Find the counterIJ for the j,k interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i)) *
                                    molality[j]*molality[k]*m_CMX_IJ_LL[counterIJ2];
                        }
                    }
                }

                // Handle neutral j species
                if (charge(j) == 0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj_LL(j,i);
                    // Zeta interaction term
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t izeta = j;
                            size_t jzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                            double zeta_LL = m_Psi_ijk_LL[n];
                            if (zeta_LL != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta_LL;
                            }
                        }
                    }
                }
            }
            // Add all of the contributions up to yield the log of the
            // solute activity coefficients (molality scale)
            m_d2lnActCoeffMolaldT2_Unscaled[i] =
                zsqd2FdT2 + sum1 + sum2 + sum3 + sum4 + sum5;
        }

        // ------ SUBSECTION FOR CALCULATING THE d2ACTCOEFFdT2 FOR ANIONS ------
        if (charge(i) < 0) {
            // species i is an anion (negative)
            double zsqd2FdT2 = charge(i)*charge(i)*d2FdT2;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                // For Anions, do the cation interactions.
                if (charge(j) > 0) {
                    sum1 += molality[j]*
                            (2.0*m_BMX_IJ_LL[counterIJ] + molarcharge*m_CMX_IJ_LL[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk_LL[n];
                            }
                        }
                    }
                }

                // For Anions, do the other anion interactions.
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ_LL[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk_LL[n];
                            // Find the counterIJ for the symmetric binary interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i)) *
                                    molality[j]*molality[k]*m_CMX_IJ_LL[counterIJ2];
                        }
                    }
                }

                // for Anions, do the neutral species interaction
                if (charge(j) == 0.0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj_LL(j,i);
                    // Zeta interaction term
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta_LL = m_Psi_ijk_LL[n];
                            if (zeta_LL != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta_LL;
                            }
                        }
                    }
                }
            }
            m_d2lnActCoeffMolaldT2_Unscaled[i] =
                zsqd2FdT2 + sum1 + sum2 + sum3 + sum4 + sum5;
        }

        // SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF
        // equations agree with my notes,
        // Equations agree with Pitzer,
        if (charge(i) == 0.0) {
            double sum1 = 0.0;
            double sum3 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                sum1 += molality[j]*2.0*m_Lambda_nj_LL(i,j);
                // Zeta term -> we piggyback on the psi term
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*m_Psi_ijk_LL[n];
                        }
                    }
                }
            }
            double sum2 = 3.0 * molality[i] * molality[i] * m_Mu_nnn_LL[i];
            m_d2lnActCoeffMolaldT2_Unscaled[i] = sum1 + sum2 + sum3;
        }
    }

    // ------ SUBSECTION FOR CALCULATING THE d2 OSMOTIC COEFF dT2 ---------
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;

    // term1 is the temperature derivative of the  DH term in the osmotic
    // coefficient expression
    // b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer implementations.
    // Is = Ionic strength on the molality scale (units of (gmol/kg))
    // Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
    double term1 = -d2AphidT2 * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (size_t j = 1; j < m_kk; j++) {
        // Loop Over Cations
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];

                    sum1 += molality[j]*molality[k] *
                            (m_BphiMX_IJ_LL[counterIJ] + molarcharge*m_CMX_IJ_LL[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between 2 cations.
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum2 += molality[j]*molality[k]*m_PhiPhi_IJ_LL[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*molality[m]*m_Psi_ijk_LL[n];
                        }
                    }
                }
            }
        }

        // Loop Over Anions
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between two anions
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];

                    sum3 += molality[j]*molality[k]*m_PhiPhi_IJ_LL[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*molality[m]*m_Psi_ijk_LL[n];
                        }
                    }
                }
            }
        }

        // Loop Over Neutral Species
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 += molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 += molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 += molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                    } else if (k == j) {
                        sum6 += 0.5 * molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            size_t n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta_LL = m_Psi_ijk_LL[n];
                            if (zeta_LL != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta_LL;
                            }
                        }
                    }
                }
            }

            sum7 += molality[j] * molality[j] * molality[j] * m_Mu_nnn_LL[j];
        }
    }
    double sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);
    // Calculate the osmotic coefficient from
    //     osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
    double d2_osmotic_coef_dT2;
    if (molalitysum > 1.0E-150) {
        d2_osmotic_coef_dT2 = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
        d2_osmotic_coef_dT2 = 0.0;
    }
    double d2_lnwateract_dT2 = -(m_weightSolvent/1000.0) * molalitysum * d2_osmotic_coef_dT2;

    // In Cantera, we define the activity coefficient of the solvent as
    //
    //     act_0 = actcoeff_0 * Xmol_0
    //
    // We have just computed act_0. However, this routine returns
    //     ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
    m_d2lnActCoeffMolaldT2_Unscaled[0] = d2_lnwateract_dT2;
}

void HMWSoln::s_update_dlnMolalityActCoeff_dP() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    m_dlnActCoeffMolaldP_Unscaled.assign(m_kk, 0.0);
    s_updatePitzer_dlnMolalityActCoeff_dP();

    for (size_t k = 1; k < m_kk; k++) {
        if (CROP_speciesCropped_[k] == 2) {
            m_dlnActCoeffMolaldP_Unscaled[k] = 0.0;
        }
    }

    if (CROP_speciesCropped_[0]) {
        m_dlnActCoeffMolaldP_Unscaled[0] = 0.0;
    }

    s_updateScaling_pHScaling_dP();
}

void HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP() const
{
    const double* molality = m_molalitiesCropped.data();

    // Local variables defined by Coltrin
    double etheta[5][5], etheta_prime[5][5], sqrtIs;

    // Molality based ionic strength of the solution
    double Is = 0.0;

    // Molar charge of the solution: In Pitzer's notation, this is his variable
    // called "Z".
    double molarcharge = 0.0;

    // molalitysum is the sum of the molalities over all solutes, even those
    // with zero charge.
    double molalitysum = 0.0;
    double currTemp = temperature();
    double currPres = pressure();

    // Make sure the counter variables are setup
    counterIJ_setup();

    // ---------- Calculate common sums over solutes ---------------------
    for (size_t n = 1; n < m_kk; n++) {
        // ionic strength
        Is += charge(n) * charge(n) * molality[n];
        // total molar charge
        molarcharge += fabs(charge(n)) * molality[n];
        molalitysum += molality[n];
    }
    Is *= 0.5;

    // Store the ionic molality in the object for reference.
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);

    // The following call to calc_lambdas() calculates all 16 elements of the
    // elambda and elambda1 arrays, given the value of the ionic strength (Is)
    calc_lambdas(Is);

    // Step 2: Find the coefficients E-theta and E-thetaprime for all
    // combinations of positive unlike charges up to 4
    for (int z1 = 1; z1 <=4; z1++) {
        for (int z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
        }
    }

    // calculate g(x) and hfunc(x) for each cation-anion pair MX
    // In the original literature, hfunc, was called gprime. However,
    // it's not the derivative of g(x), so I renamed it.
    for (size_t i = 1; i < (m_kk - 1); i++) {
        for (size_t j = (i+1); j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // Only loop over oppositely charge species
            if (charge(i)*charge(j) < 0) {
                // x is a reduced function variable
                double x1 = sqrtIs * m_Alpha1MX_ij[counterIJ];
                if (x1 > 1.0E-100) {
                    m_gfunc_IJ[counterIJ] = 2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 * x1);
                    m_hfunc_IJ[counterIJ] = -2.0*
                                       (1.0-(1.0 + x1 + 0.5 * x1 * x1) * exp(-x1)) / (x1 * x1);
                } else {
                    m_gfunc_IJ[counterIJ] = 0.0;
                    m_hfunc_IJ[counterIJ] = 0.0;
                }

                if (m_Beta2MX_ij_P[counterIJ] != 0.0) {
                    double x2 = sqrtIs * m_Alpha2MX_ij[counterIJ];
                    if (x2 > 1.0E-100) {
                        m_g2func_IJ[counterIJ] = 2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        m_h2func_IJ[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        m_g2func_IJ[counterIJ] = 0.0;
                        m_h2func_IJ[counterIJ] = 0.0;
                    }
                }
            } else {
                m_gfunc_IJ[counterIJ] = 0.0;
                m_hfunc_IJ[counterIJ] = 0.0;
            }
        }
    }

    // SUBSECTION TO CALCULATE BMX_P, BprimeMX_P, BphiMX_P
    // These are now temperature derivatives of the previously calculated
    // quantities.
    for (size_t i = 1; i < m_kk - 1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_BMX_IJ_P[counterIJ] = m_Beta0MX_ij_P[counterIJ]
                                    + m_Beta1MX_ij_P[counterIJ] * m_gfunc_IJ[counterIJ]
                                    + m_Beta2MX_ij_P[counterIJ] * m_g2func_IJ[counterIJ];
                if (Is > 1.0E-150) {
                    m_BprimeMX_IJ_P[counterIJ] = (m_Beta1MX_ij_P[counterIJ] * m_hfunc_IJ[counterIJ]/Is +
                                             m_Beta2MX_ij_P[counterIJ] * m_h2func_IJ[counterIJ]/Is);
                } else {
                    m_BprimeMX_IJ_P[counterIJ] = 0.0;
                }
                m_BphiMX_IJ_P[counterIJ] = m_BMX_IJ_P[counterIJ] + Is*m_BprimeMX_IJ_P[counterIJ];
            } else {
                m_BMX_IJ_P[counterIJ] = 0.0;
                m_BprimeMX_IJ_P[counterIJ] = 0.0;
                m_BphiMX_IJ_P[counterIJ] = 0.0;
            }
        }
    }

    // --------- SUBSECTION TO CALCULATE CMX_P ----------
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0.0) {
                m_CMX_IJ_P[counterIJ] = m_CphiMX_ij_P[counterIJ]/
                                   (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                m_CMX_IJ_P[counterIJ] = 0.0;
            }
        }
    }

    // ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) > 0) {
                m_Phi_IJ_P[counterIJ] = m_Theta_ij_P[counterIJ];
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ_P[counterIJ] = m_Phi_IJ_P[counterIJ] + Is * m_Phiprime_IJ[counterIJ];
            } else {
                m_Phi_IJ_P[counterIJ] = 0.0;
                m_Phiprime_IJ[counterIJ] = 0.0;
                m_PhiPhi_IJ_P[counterIJ] = 0.0;
            }
        }
    }

    // ----------- SUBSECTION FOR CALCULATION OF dFdT ---------------------
    double dA_DebyedP = dA_DebyedP_TP(currTemp, currPres);
    double dAphidP = dA_DebyedP /3.0;
    double dFdP = -dAphidP * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                       + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    for (size_t i = 1; i < m_kk-1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            // Find the counterIJ for the symmetric binary interaction
            size_t n = m_kk*i + j;
            size_t counterIJ = m_CounterIJ[n];

            // both species have a non-zero charge, and one is positive
            // and the other is negative
            if (charge(i)*charge(j) < 0) {
                dFdP += molality[i]*molality[j] * m_BprimeMX_IJ_P[counterIJ];
            }

            // Both species have a non-zero charge, and they
            // have the same sign, e.g., both positive or both negative.
            if (charge(i)*charge(j) > 0) {
                dFdP += molality[i]*molality[j] * m_Phiprime_IJ[counterIJ];
            }
        }
    }

    for (size_t i = 1; i < m_kk; i++) {
        // -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdP FOR CATIONS -----
        if (charge(i) > 0) {
            // species i is the cation (positive) to calc the actcoeff
            double zsqdFdP = charge(i)*charge(i)*dFdP;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 += molality[j]*
                            (2.0*m_BMX_IJ_P[counterIJ] + molarcharge*m_CMX_IJ_P[counterIJ]);
                    if (j < m_kk-1) {
                        // This term is the ternary interaction involving the
                        // non-duplicate sum over double anions, j, k, with
                        // respect to the cation, i.
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk_P[n];
                            }
                        }
                    }
                }

                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ_P[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk_P[n];

                            // Find the counterIJ for the j,k interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i)) *
                                    molality[j]*molality[k]*m_CMX_IJ_P[counterIJ2];
                        }
                    }
                }

                // for Anions, do the neutral species interaction
                if (charge(j) == 0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj_P(j,i);
                    // Zeta interaction term
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t izeta = j;
                            size_t jzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                            double zeta_P = m_Psi_ijk_P[n];
                            if (zeta_P != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta_P;
                            }
                        }
                    }
                }
            }

            // Add all of the contributions up to yield the log of the
            // solute activity coefficients (molality scale)
            m_dlnActCoeffMolaldP_Unscaled[i] =
                zsqdFdP + sum1 + sum2 + sum3 + sum4 + sum5;
        }

        // ------ SUBSECTION FOR CALCULATING THE dACTCOEFFdP FOR ANIONS ------
        if (charge(i) < 0) {
            // species i is an anion (negative)
            double zsqdFdP = charge(i)*charge(i)*dFdP;
            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                // Find the counterIJ for the symmetric binary interaction
                size_t n = m_kk*i + j;
                size_t counterIJ = m_CounterIJ[n];

                // For Anions, do the cation interactions.
                if (charge(j) > 0) {
                    sum1 += molality[j] *
                            (2.0*m_BMX_IJ_P[counterIJ] + molarcharge*m_CMX_IJ_P[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 += molality[j]*molality[k]*m_Psi_ijk_P[n];
                            }
                        }
                    }
                }

                // For Anions, do the other anion interactions.
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 += molality[j]*(2.0*m_Phi_IJ_P[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*m_Psi_ijk_P[n];
                            // Find the counterIJ for the symmetric binary interaction
                            n = m_kk*j + k;
                            size_t counterIJ2 = m_CounterIJ[n];
                            sum4 += fabs(charge(i))*
                                    molality[j]*molality[k]*m_CMX_IJ_P[counterIJ2];
                        }
                    }
                }

                // for Anions, do the neutral species interaction
                if (charge(j) == 0.0) {
                    sum5 += molality[j]*2.0*m_Lambda_nj_P(j,i);
                    // Zeta interaction term
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta_P = m_Psi_ijk_P[n];
                            if (zeta_P != 0.0) {
                                sum5 += molality[j]*molality[k]*zeta_P;
                            }
                        }
                    }
                }
            }
            m_dlnActCoeffMolaldP_Unscaled[i] =
                zsqdFdP + sum1 + sum2 + sum3 + sum4 + sum5;
        }

        // ------ SUBSECTION FOR CALCULATING d NEUTRAL SOLUTE ACT COEFF dP -----
        if (charge(i) == 0.0) {
            double sum1 = 0.0;
            double sum3 = 0.0;
            for (size_t j = 1; j < m_kk; j++) {
                sum1 += molality[j]*2.0*m_Lambda_nj_P(i,j);
                // Zeta term -> we piggyback on the psi term
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*m_Psi_ijk_P[n];
                        }
                    }
                }
            }
            double sum2 = 3.0 * molality[i] * molality[i] * m_Mu_nnn_P[i];
            m_dlnActCoeffMolaldP_Unscaled[i] = sum1 + sum2 + sum3;
        }
    }

    // ------ SUBSECTION FOR CALCULATING THE d OSMOTIC COEFF dP ---------
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;

    // term1 is the temperature derivative of the DH term in the osmotic
    // coefficient expression
    // b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer implementations.
    // Is = Ionic strength on the molality scale (units of (gmol/kg))
    // Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
    double term1 = -dAphidP * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (size_t j = 1; j < m_kk; j++) {
        // Loop Over Cations
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum1 += molality[j]*molality[k]*
                            (m_BphiMX_IJ_P[counterIJ] + molarcharge*m_CMX_IJ_P[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between 2 cations.
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];
                    sum2 += molality[j]*molality[k]*m_PhiPhi_IJ_P[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 += molality[j]*molality[k]*molality[m]*m_Psi_ijk_P[n];
                        }
                    }
                }
            }
        }

        // Loop Over Anions
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    // Find the counterIJ for the symmetric j,k binary interaction
                    // between two anions
                    size_t n = m_kk*j + k;
                    size_t counterIJ = m_CounterIJ[n];

                    sum3 += molality[j]*molality[k]*m_PhiPhi_IJ_P[counterIJ];
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 += molality[j]*molality[k]*molality[m]*m_Psi_ijk_P[n];
                        }
                    }
                }
            }
        }

        // Loop Over Neutral Species
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 += molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 += molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 += molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                    } else if (k == j) {
                        sum6 += 0.5 * molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (size_t m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            size_t n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta_P = m_Psi_ijk_P[n];
                            if (zeta_P != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta_P;
                            }
                        }
                    }
                }
            }

            sum7 += molality[j] * molality[j] * molality[j] * m_Mu_nnn_P[j];
        }
    }
    double sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);

    // Calculate the osmotic coefficient from
    //     osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
    double d_osmotic_coef_dP;
    if (molalitysum > 1.0E-150) {
        d_osmotic_coef_dP = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
        d_osmotic_coef_dP = 0.0;
    }
    double d_lnwateract_dP = -(m_weightSolvent/1000.0) * molalitysum * d_osmotic_coef_dP;

    // In Cantera, we define the activity coefficient of the solvent as
    //
    //     act_0 = actcoeff_0 * Xmol_0
    //
    // We have just computed act_0. However, this routine returns
    //     ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
    m_dlnActCoeffMolaldP_Unscaled[0] = d_lnwateract_dP;
}

void HMWSoln::calc_lambdas(double is) const
{
    if( m_last_is == is ) {
      return;
    }
    m_last_is = is;

    // Coefficients c1-c4 are used to approximate the integral function "J";
    // aphi is the Debye-Huckel constant at 25 C
    double c1 = 4.581, c2 = 0.7237, c3 = 0.0120, c4 = 0.528;
    double aphi = 0.392; /* Value at 25 C */
    if (is < 1.0E-150) {
        for (int i = 0; i < 17; i++) {
            elambda[i] = 0.0;
            elambda1[i] = 0.0;
        }
        return;
    }

    // Calculate E-lambda terms for charge combinations of like sign,
    // using method of Pitzer (1975). Charges up to 4 are calculated.
    for (int i=1; i<=4; i++) {
        for (int j=i; j<=4; j++) {
            int ij = i*j;

            // calculate the product of the charges
            double zprod = (double)ij;

            // calculate Xmn (A1) from Harvie, Weare (1980).
            double x = 6.0* zprod * aphi * sqrt(is); // eqn 23

            double jfunc = x / (4.0 + c1*pow(x,-c2)*exp(-c3*pow(x,c4))); // eqn 47

            double t = c3 * c4 * pow(x,c4);
            double dj = c1* pow(x,(-c2-1.0)) * (c2+t) * exp(-c3*pow(x,c4));
            double jprime = (jfunc/x)*(1.0 + jfunc*dj);

            elambda[ij] = zprod*jfunc / (4.0*is); // eqn 14
            elambda1[ij] = (3.0*zprod*zprod*aphi*jprime/(4.0*sqrt(is))
                            - elambda[ij])/is;
        }
    }
}

void HMWSoln::calc_thetas(int z1, int z2,
                          double* etheta, double* etheta_prime) const
{
    // Calculate E-theta(i) and E-theta'(I) using method of Pitzer (1987)
    int i = abs(z1);
    int j = abs(z2);

    AssertThrowMsg(i <= 4 && j <= 4, "HMWSoln::calc_thetas",
                   "we shouldn't be here");
    AssertThrowMsg(i != 0 && j != 0, "HMWSoln::calc_thetas",
                   "called with one species being neutral");

    // Check to see if the charges are of opposite sign. If they are of opposite
    // sign then their etheta interaction is zero.
    if (z1*z2 < 0) {
        *etheta = 0.0;
        *etheta_prime = 0.0;
    } else {
        // Actually calculate the interaction.
        double f1 = (double)i / (2.0 * j);
        double f2 = (double)j / (2.0 * i);
        *etheta = elambda[i*j] - f1*elambda[j*j] - f2*elambda[i*i];
        *etheta_prime = elambda1[i*j] - f1*elambda1[j*j] - f2*elambda1[i*i];
    }
}

void HMWSoln::s_updateIMS_lnMolalityActCoeff() const
{
    // Calculate the molalities. Currently, the molalities may not be current
    // with respect to the contents of the State objects' data.
    calcMolalities();
    double xmolSolvent = moleFraction(0);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);
    // Exponentials - trial 2
    if (xmolSolvent > IMS_X_o_cutoff_) {
        for (size_t k = 1; k < m_kk; k++) {
            IMS_lnActCoeffMolal_[k]= 0.0;
        }
        IMS_lnActCoeffMolal_[0] = - log(xx) + (xx - 1.0)/xx;
        return;
    } else {
        double xoverc = xmolSolvent/IMS_cCut_;
        double eterm = std::exp(-xoverc);

        double fptmp = IMS_bfCut_ - IMS_afCut_ / IMS_cCut_ - IMS_bfCut_*xoverc
                       + 2.0*IMS_dfCut_*xmolSolvent - IMS_dfCut_*xmolSolvent*xoverc;
        double f_prime = 1.0 + eterm*fptmp;
        double f = xmolSolvent + IMS_efCut_
                   + eterm * (IMS_afCut_ + xmolSolvent * (IMS_bfCut_ + IMS_dfCut_*xmolSolvent));

        double gptmp = IMS_bgCut_ - IMS_agCut_ / IMS_cCut_ - IMS_bgCut_*xoverc
                       + 2.0*IMS_dgCut_*xmolSolvent - IMS_dgCut_*xmolSolvent*xoverc;
        double g_prime = 1.0 + eterm*gptmp;
        double g = xmolSolvent + IMS_egCut_
                   + eterm * (IMS_agCut_ + xmolSolvent * (IMS_bgCut_ + IMS_dgCut_*xmolSolvent));

        double tmp = (xmolSolvent / g * g_prime + (1.0 - xmolSolvent) / f * f_prime);
        double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
        double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

        tmp = log(xx) + lngammak;
        for (size_t k = 1; k < m_kk; k++) {
            IMS_lnActCoeffMolal_[k]= tmp;
        }
        IMS_lnActCoeffMolal_[0] = lngammao;
    }
    return;
}

void HMWSoln::printCoeffs() const
{
    calcMolalities();
    vector_fp& moleF = m_tmpV;

    // Update the coefficients wrt Temperature. Calculate the derivatives as well
    s_updatePitzer_CoeffWRTemp(2);
    getMoleFractions(moleF.data());

    writelog("Index  Name                  MoleF   MolalityCropped  Charge\n");
    for (size_t k = 0; k < m_kk; k++) {
        writelogf("%2d     %-16s %14.7le %14.7le %5.1f \n",
                  k, speciesName(k), moleF[k], m_molalitiesCropped[k], charge(k));
    }

    writelog("\n Species          Species            beta0MX  "
           "beta1MX   beta2MX   CphiMX    alphaMX thetaij\n");
    for (size_t i = 1; i < m_kk - 1; i++) {
        for (size_t j = i+1; j < m_kk; j++) {
            size_t n = i * m_kk + j;
            size_t ct = m_CounterIJ[n];
            writelogf(" %-16s %-16s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n",
                      speciesName(i), speciesName(j),
                      m_Beta0MX_ij[ct], m_Beta1MX_ij[ct],
                      m_Beta2MX_ij[ct], m_CphiMX_ij[ct],
                      m_Alpha1MX_ij[ct], m_Theta_ij[ct]);
        }
    }

    writelog("\n Species          Species          Species       psi   \n");
    for (size_t i = 1; i < m_kk; i++) {
        for (size_t j = 1; j < m_kk; j++) {
            for (size_t k = 1; k < m_kk; k++) {
                size_t n = k + j * m_kk + i * m_kk * m_kk;
                if (m_Psi_ijk[n] != 0.0) {
                    writelogf(" %-16s %-16s %-16s %9.5f \n",
                              speciesName(i), speciesName(j),
                              speciesName(k), m_Psi_ijk[n]);
                }
            }
        }
    }
}

void HMWSoln::applyphScale(doublereal* acMolality) const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal lnGammaClMs2 = s_NBS_CLM_lnMolalityActCoeff();
    doublereal lnGammaCLMs1 = m_lnActCoeffMolal_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(lnGammaClMs2 - lnGammaCLMs1);
    for (size_t k = 0; k < m_kk; k++) {
        acMolality[k] *= exp(charge(k) * afac);
    }
}

void HMWSoln::s_updateScaling_pHScaling() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_lnActCoeffMolal_Scaled = m_lnActCoeffMolal_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal lnGammaClMs2 = s_NBS_CLM_lnMolalityActCoeff();
    doublereal lnGammaCLMs1 = m_lnActCoeffMolal_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(lnGammaClMs2 - lnGammaCLMs1);
    for (size_t k = 0; k < m_kk; k++) {
        m_lnActCoeffMolal_Scaled[k] = m_lnActCoeffMolal_Unscaled[k] + charge(k) * afac;
    }
}

void HMWSoln::s_updateScaling_pHScaling_dT() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_dlnActCoeffMolaldT_Scaled = m_dlnActCoeffMolaldT_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal dlnGammaClM_dT_s2 = s_NBS_CLM_dlnMolalityActCoeff_dT();
    doublereal dlnGammaCLM_dT_s1 = m_dlnActCoeffMolaldT_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(dlnGammaClM_dT_s2 - dlnGammaCLM_dT_s1);
    for (size_t k = 0; k < m_kk; k++) {
        m_dlnActCoeffMolaldT_Scaled[k] = m_dlnActCoeffMolaldT_Unscaled[k] + charge(k) * afac;
    }
}

void HMWSoln::s_updateScaling_pHScaling_dT2() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_d2lnActCoeffMolaldT2_Scaled = m_d2lnActCoeffMolaldT2_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal d2lnGammaClM_dT2_s2 = s_NBS_CLM_d2lnMolalityActCoeff_dT2();
    doublereal d2lnGammaCLM_dT2_s1 = m_d2lnActCoeffMolaldT2_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(d2lnGammaClM_dT2_s2 - d2lnGammaCLM_dT2_s1);
    for (size_t k = 0; k < m_kk; k++) {
        m_d2lnActCoeffMolaldT2_Scaled[k] = m_d2lnActCoeffMolaldT2_Unscaled[k] + charge(k) * afac;
    }
}

void HMWSoln::s_updateScaling_pHScaling_dP() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_dlnActCoeffMolaldP_Scaled = m_dlnActCoeffMolaldP_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal dlnGammaClM_dP_s2 = s_NBS_CLM_dlnMolalityActCoeff_dP();
    doublereal dlnGammaCLM_dP_s1 = m_dlnActCoeffMolaldP_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(dlnGammaClM_dP_s2 - dlnGammaCLM_dP_s1);
    for (size_t k = 0; k < m_kk; k++) {
        m_dlnActCoeffMolaldP_Scaled[k] = m_dlnActCoeffMolaldP_Unscaled[k] + charge(k) * afac;
    }
}

doublereal HMWSoln::s_NBS_CLM_lnMolalityActCoeff() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal A = A_Debye_TP();
    doublereal lnGammaClMs2 = - A * sqrtIs /(1.0 + 1.5 * sqrtIs);
    return lnGammaClMs2;
}

doublereal HMWSoln::s_NBS_CLM_dlnMolalityActCoeff_dT() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal dAdT = dA_DebyedT_TP();
    return - dAdT * sqrtIs /(1.0 + 1.5 * sqrtIs);
}

doublereal HMWSoln::s_NBS_CLM_d2lnMolalityActCoeff_dT2() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal d2AdT2 = d2A_DebyedT2_TP();
    return - d2AdT2 * sqrtIs /(1.0 + 1.5 * sqrtIs);
}

doublereal HMWSoln::s_NBS_CLM_dlnMolalityActCoeff_dP() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal dAdP = dA_DebyedP_TP();
    return - dAdP * sqrtIs /(1.0 + 1.5 * sqrtIs);
}

}
