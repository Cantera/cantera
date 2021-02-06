/**
 *  @file IdealMolalSoln.cpp
 *   ThermoPhase object for the ideal molal equation of
 * state (see \ref thermoprops
 * and class \link Cantera::IdealMolalSoln IdealMolalSoln\endlink).
 *
 * Definition file for a derived class of ThermoPhase that handles variable
 * pressure standard state methods for calculating thermodynamic properties that
 * are further based upon activities on the molality scale. The Ideal molal
 * solution assumes that all molality-based activity coefficients are equal to
 * one. This turns out, actually, to be highly nonlinear when the solvent
 * densities get low.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

#include <iostream>

namespace Cantera
{

IdealMolalSoln::IdealMolalSoln() :
    m_formGC(2),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.20),
    IMS_gamma_o_min_(0.00001),
    IMS_gamma_k_min_(10.0),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_cCut_(.05),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0)
{
}

IdealMolalSoln::IdealMolalSoln(const std::string& inputFile,
                               const std::string& id_) :
    MolalityVPSSTP(),
    m_formGC(2),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(0.00001),
    IMS_gamma_k_min_(10.0),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_cCut_(.05),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0)
{
    initThermoFile(inputFile, id_);
}

IdealMolalSoln::IdealMolalSoln(XML_Node& root, const std::string& id_) :
    MolalityVPSSTP(),
    m_formGC(2),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(0.00001),
    IMS_gamma_k_min_(10.0),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_cCut_(.05),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0)
{
    importPhase(root, this);
}

doublereal IdealMolalSoln::enthalpy_mole() const
{
    getPartialMolarEnthalpies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal IdealMolalSoln::intEnergy_mole() const
{
    getPartialMolarEnthalpies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal IdealMolalSoln::entropy_mole() const
{
    getPartialMolarEntropies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal IdealMolalSoln::gibbs_mole() const
{
    getChemPotentials(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal IdealMolalSoln::cp_mole() const
{
    getPartialMolarCp(m_tmpV.data());
    return mean_X(m_tmpV);
}

// ------- Mechanical Equation of State Properties ------------------------

void IdealMolalSoln::calcDensity()
{
    getPartialMolarVolumes(m_tmpV.data());
    doublereal dd = meanMolecularWeight() / mean_X(m_tmpV);
    Phase::assignDensity(dd);
}

doublereal IdealMolalSoln::isothermalCompressibility() const
{
    return 0.0;
}

doublereal IdealMolalSoln::thermalExpansionCoeff() const
{
    return 0.0;
}

void IdealMolalSoln::setDensity(const doublereal rho)
{
    warn_deprecated("IdealMolalSoln::setDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setDensity instead");
    if (rho != density()) {
        throw CanteraError("IdealMolalSoln::setDensity",
                           "Density is not an independent variable");
    }
}

void IdealMolalSoln::setMolarDensity(const doublereal conc)
{
    warn_deprecated("IdealMolalSoln::setMolarDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setMolarDensity instead");
    if (conc != Phase::molarDensity()) {
        throw CanteraError("IdealMolalSoln::setMolarDensity",
                           "molarDensity/density is not an independent variable");
    }
}

// ------- Activities and Activity Concentrations

Units IdealMolalSoln::standardConcentrationUnits() const
{
    if (m_formGC == 0) {
        return Units(1.0); // dimensionless
    } else {
        // kmol/m^3 for bulk phases
        return Units(1.0, 0, -static_cast<double>(nDim()), 0, 0, 0, 1);
    }
}

void IdealMolalSoln::getActivityConcentrations(doublereal* c) const
{
    if (m_formGC != 1) {
        double c_solvent = standardConcentration();
        getActivities(c);
        for (size_t k = 0; k < m_kk; k++) {
            c[k] *= c_solvent;
        }
    } else {
        getActivities(c);
        for (size_t k = 0; k < m_kk; k++) {
            double c0 = standardConcentration(k);
            c[k] *= c0;
        }
    }
}

doublereal IdealMolalSoln::standardConcentration(size_t k) const
{
    switch (m_formGC) {
    case 0:
        return 1.0;
    case 1:
        return 1.0 / m_speciesMolarVolume[k];
    case 2:
        return 1.0 / m_speciesMolarVolume[0];
    default:
        throw CanteraError("IdealMolalSoln::standardConcentration",
                       "m_formGC is set to an incorrect value. \
                       Allowed values are 0, 1, and 2");
    }
}

void IdealMolalSoln::getActivities(doublereal* ac) const
{
    _updateStandardStateThermo();

    // Update the molality array, m_molalities(). This requires an update due to
    // mole fractions
    if (IMS_typeCutoff_ == 0) {
        calcMolalities();
        for (size_t k = 0; k < m_kk; k++) {
            ac[k] = m_molalities[k];
        }
        double xmolSolvent = moleFraction(0);
        // Limit the activity coefficient to be finite as the solvent mole
        // fraction goes to zero.
        xmolSolvent = std::max(m_xmolSolventMIN, xmolSolvent);
        ac[0] = exp((xmolSolvent - 1.0)/xmolSolvent);
    } else {

        s_updateIMS_lnMolalityActCoeff();

        // Now calculate the array of activities.
        for (size_t k = 1; k < m_kk; k++) {
            ac[k] = m_molalities[k] * exp(IMS_lnActCoeffMolal_[k]);
        }
        double xmolSolvent = moleFraction(0);
        ac[0] = exp(IMS_lnActCoeffMolal_[0]) * xmolSolvent;
    }
}

void IdealMolalSoln::getMolalityActivityCoefficients(doublereal* acMolality) const
{
    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 0; k < m_kk; k++) {
            acMolality[k] = 1.0;
        }
        double xmolSolvent = moleFraction(0);
        // Limit the activity coefficient to be finite as the solvent mole
        // fraction goes to zero.
        xmolSolvent = std::max(m_xmolSolventMIN, xmolSolvent);
        acMolality[0] = exp((xmolSolvent - 1.0)/xmolSolvent) / xmolSolvent;
    } else {
        s_updateIMS_lnMolalityActCoeff();
        std::copy(IMS_lnActCoeffMolal_.begin(), IMS_lnActCoeffMolal_.end(), acMolality);
        for (size_t k = 0; k < m_kk; k++) {
            acMolality[k] = exp(acMolality[k]);
        }
    }
}

// ------ Partial Molar Properties of the Solution -----------------

void IdealMolalSoln::getChemPotentials(doublereal* mu) const
{
    // First get the standard chemical potentials. This requires updates of
    // standard state as a function of T and P These are defined at unit
    // molality.
    getStandardChemPotentials(mu);

    // Update the molality array, m_molalities(). This requires an update due to
    // mole fractions
    calcMolalities();

    // get the solvent mole fraction
    double xmolSolvent = moleFraction(0);

    if (IMS_typeCutoff_ == 0 || xmolSolvent > 3.* IMS_X_o_cutoff_/2.0) {
        for (size_t k = 1; k < m_kk; k++) {
            double xx = std::max(m_molalities[k], SmallNumber);
            mu[k] += RT() * log(xx);
        }

        // Do the solvent
        //  -> see my notes
        double xx = std::max(xmolSolvent, SmallNumber);
        mu[0] += (RT() * (xmolSolvent - 1.0) / xx);
    } else {
        // Update the activity coefficients. This also updates the internal
        // molality array.
        s_updateIMS_lnMolalityActCoeff();

        for (size_t k = 1; k < m_kk; k++) {
            double xx = std::max(m_molalities[k], SmallNumber);
            mu[k] += RT() * (log(xx) + IMS_lnActCoeffMolal_[k]);
        }
        double xx = std::max(xmolSolvent, SmallNumber);
        mu[0] += RT() * (log(xx) + IMS_lnActCoeffMolal_[0]);
    }
}

void IdealMolalSoln::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }
}

void IdealMolalSoln::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    calcMolalities();
    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 1; k < m_kk; k++) {
            doublereal mm = std::max(SmallNumber, m_molalities[k]);
            sbar[k] -= GasConstant * log(mm);
        }
        double xmolSolvent = moleFraction(0);
        sbar[0] -= (GasConstant * (xmolSolvent - 1.0) / xmolSolvent);
    } else {
        // Update the activity coefficients, This also update the internally
        // stored molalities.
        s_updateIMS_lnMolalityActCoeff();

        // First we will add in the obvious dependence on the T term out front
        // of the log activity term
        doublereal mm;
        for (size_t k = 1; k < m_kk; k++) {
            mm = std::max(SmallNumber, m_molalities[k]);
            sbar[k] -= GasConstant * (log(mm) + IMS_lnActCoeffMolal_[k]);
        }
        double xmolSolvent = moleFraction(0);
        mm = std::max(SmallNumber, xmolSolvent);
        sbar[0] -= GasConstant *(log(mm) + IMS_lnActCoeffMolal_[0]);
    }
}

void IdealMolalSoln::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

void IdealMolalSoln::getPartialMolarCp(doublereal* cpbar) const
{
    // Get the nondimensional Gibbs standard state of the species at the T and P
    // of the solution.
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

// -------------- Utilities -------------------------------

bool IdealMolalSoln::addSpecies(shared_ptr<Species> spec)
{
    bool added = MolalityVPSSTP::addSpecies(spec);
    if (added) {
        m_speciesMolarVolume.push_back(0.0);
        m_tmpV.push_back(0.0);
        IMS_lnActCoeffMolal_.push_back(0.0);
    }
    return added;
}

void IdealMolalSoln::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    MolalityVPSSTP::initThermoXML(phaseNode, id_);

    if (id_.size() > 0 && phaseNode.id() != id_) {
        throw CanteraError("IdealMolalSoln::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("IdealMolalSoln::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    // Possible change the form of the standard concentrations
    if (thermoNode.hasChild("standardConc")) {
        XML_Node& scNode = thermoNode.child("standardConc");
        setStandardConcentrationModel(scNode["model"]);
    }

    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        std::string modelString = acNode.attrib("model");
        if (modelString != "IdealMolalSoln") {
            throw CanteraError("IdealMolalSoln::initThermoXML",
                               "unknown ActivityCoefficient model: " + modelString);
        }
        if (acNode.hasChild("idealMolalSolnCutoff")) {
            XML_Node& ccNode = acNode.child("idealMolalSolnCutoff");
            modelString = ccNode.attrib("model");
            if (modelString != "") {
                setCutoffModel(modelString);
                if (ccNode.hasChild("gamma_o_limit")) {
                    IMS_gamma_o_min_ = getFloat(ccNode, "gamma_o_limit");
                }
                if (ccNode.hasChild("gamma_k_limit")) {
                    IMS_gamma_k_min_ = getFloat(ccNode, "gamma_k_limit");
                }
                if (ccNode.hasChild("X_o_cutoff")) {
                    IMS_X_o_cutoff_ = getFloat(ccNode, "X_o_cutoff");
                }
                if (ccNode.hasChild("c_0_param")) {
                    IMS_cCut_ = getFloat(ccNode, "c_0_param");
                }
                if (ccNode.hasChild("slope_f_limit")) {
                    IMS_slopefCut_ = getFloat(ccNode, "slope_f_limit");
                }
                if (ccNode.hasChild("slope_g_limit")) {
                    IMS_slopegCut_ = getFloat(ccNode, "slope_g_limit");
                }
            }
        } else {
            setCutoffModel("none");
        }
    }
}

void IdealMolalSoln::initThermo()
{
    MolalityVPSSTP::initThermo();

    if (m_input.hasKey("standard-concentration-basis")) {
        setStandardConcentrationModel(m_input["standard-concentration-basis"].asString());
    }
    if (m_input.hasKey("cutoff")) {
        auto& cutoff = m_input["cutoff"].as<AnyMap>();
        setCutoffModel(cutoff.getString("model", "none"));
        if (cutoff.hasKey("gamma_o")) {
            IMS_gamma_o_min_ = cutoff["gamma_o"].asDouble();
        }
        if (cutoff.hasKey("gamma_k")) {
            IMS_gamma_k_min_ = cutoff["gamma_k"].asDouble();
        }
        if (cutoff.hasKey("X_o")) {
            IMS_X_o_cutoff_ = cutoff["X_o"].asDouble();
        }
        if (cutoff.hasKey("c_0")) {
            IMS_cCut_ = cutoff["c_0"].asDouble();
        }
        if (cutoff.hasKey("slope_f")) {
            IMS_slopefCut_ = cutoff["slope_f"].asDouble();
        }
        if (cutoff.hasKey("slope_g")) {
            IMS_slopegCut_ = cutoff["slope_g"].asDouble();
        }
    }

    for (size_t k = 0; k < nSpecies(); k++) {
        m_speciesMolarVolume[k] = providePDSS(k)->molarVolume();
    }
    if (IMS_typeCutoff_ == 2) {
        calcIMSCutoffParams_();
    }
    setMoleFSolventMin(1.0E-5);
}

void IdealMolalSoln::setStandardConcentrationModel(const std::string& model)
{
    if (caseInsensitiveEquals(model, "unity")) {
        m_formGC = 0;
    } else if (caseInsensitiveEquals(model, "species-molar-volume")
               || caseInsensitiveEquals(model, "molar_volume")) {
        m_formGC = 1;
    } else if (caseInsensitiveEquals(model, "solvent-molar-volume")
               || caseInsensitiveEquals(model, "solvent_volume")) {
        m_formGC = 2;
    } else {
        throw CanteraError("IdealMolalSoln::setStandardConcentrationModel",
                           "Unknown standard concentration model '{}'", model);
    }
}

void IdealMolalSoln::setCutoffModel(const std::string& model)
{
    if (caseInsensitiveEquals(model, "none")) {
        IMS_typeCutoff_ = 0;
    } else if (caseInsensitiveEquals(model, "poly")) {
        IMS_typeCutoff_ = 1;
    } else if (caseInsensitiveEquals(model, "polyexp")) {
        IMS_typeCutoff_ = 2;
    } else {
        throw CanteraError("IdealMolalSoln::setCutoffModel",
                           "Unknown cutoff model '{}'", model);
    }
}

// ------------ Private and Restricted Functions ------------------

void IdealMolalSoln::s_updateIMS_lnMolalityActCoeff() const
{
    // Calculate the molalities. Currently, the molalities may not be current
    // with respect to the contents of the State objects' data.
    calcMolalities();

    double xmolSolvent = moleFraction(0);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);

    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 1; k < m_kk; k++) {
            IMS_lnActCoeffMolal_[k]= 0.0;
        }
        IMS_lnActCoeffMolal_[0] = - log(xx) + (xx - 1.0)/xx;
        return;
    } else if (IMS_typeCutoff_ == 1) {
        if (xmolSolvent > 3.0 * IMS_X_o_cutoff_/2.0) {
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= 0.0;
            }
            IMS_lnActCoeffMolal_[0] = - log(xx) + (xx - 1.0)/xx;
            return;
        } else if (xmolSolvent < IMS_X_o_cutoff_/2.0) {
            double tmp = log(xx * IMS_gamma_k_min_);
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[0] = log(IMS_gamma_o_min_);
            return;
        } else {
            // If we are in the middle region, calculate the connecting polynomials
            double xminus = xmolSolvent - IMS_X_o_cutoff_/2.0;
            double xminus2 = xminus * xminus;
            double xminus3 = xminus2 * xminus;
            double x_o_cut2 = IMS_X_o_cutoff_ * IMS_X_o_cutoff_;
            double x_o_cut3 = x_o_cut2 * IMS_X_o_cutoff_;

            double h2 = 3.5 * xminus2 / IMS_X_o_cutoff_ - 2.0 * xminus3 / x_o_cut2;
            double h2_prime = 7.0 * xminus / IMS_X_o_cutoff_ - 6.0 * xminus2 / x_o_cut2;

            double h1 = (1.0 - 3.0 * xminus2 / x_o_cut2 + 2.0 * xminus3/ x_o_cut3);
            double h1_prime = (- 6.0 * xminus / x_o_cut2 + 6.0 * xminus2/ x_o_cut3);

            double h1_g = h1 / IMS_gamma_o_min_;
            double h1_g_prime = h1_prime / IMS_gamma_o_min_;

            double alpha = 1.0 / (exp(1.0) * IMS_gamma_k_min_);
            double h1_f = h1 * alpha;
            double h1_f_prime = h1_prime * alpha;

            double f = h2 + h1_f;
            double f_prime = h2_prime + h1_f_prime;

            double g = h2 + h1_g;
            double g_prime = h2_prime + h1_g_prime;

            double tmp = (xmolSolvent/ g * g_prime + (1.0-xmolSolvent) / f * f_prime);
            double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
            double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

            tmp = log(xmolSolvent) + lngammak;
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[0] = lngammao;
        }
    } else if (IMS_typeCutoff_ == 2) {
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
            double f = xmolSolvent + IMS_efCut_ + eterm * (IMS_afCut_ + xmolSolvent * (IMS_bfCut_ + IMS_dfCut_*xmolSolvent));

            double gptmp = IMS_bgCut_ - IMS_agCut_ / IMS_cCut_ - IMS_bgCut_*xoverc
                           + 2.0*IMS_dgCut_*xmolSolvent - IMS_dgCut_*xmolSolvent*xoverc;
            double g_prime = 1.0 + eterm*gptmp;
            double g = xmolSolvent + IMS_egCut_ + eterm * (IMS_agCut_ + xmolSolvent * (IMS_bgCut_ + IMS_dgCut_*xmolSolvent));

            double tmp = (xmolSolvent / g * g_prime + (1.0 - xmolSolvent) / f * f_prime);
            double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
            double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

            tmp = log(xx) + lngammak;
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[0] = lngammao;
        }
    }
}

void IdealMolalSoln::calcIMSCutoffParams_()
{
    IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_);
    IMS_efCut_ = 0.0;
    bool converged = false;
    for (int its = 0; its < 100 && !converged; its++) {
        double oldV = IMS_efCut_;
        IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_) - IMS_efCut_;
        IMS_bfCut_ = IMS_afCut_ / IMS_cCut_ + IMS_slopefCut_ - 1.0;
        IMS_dfCut_ = ((- IMS_afCut_/IMS_cCut_ + IMS_bfCut_ - IMS_bfCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_afCut_ + IMS_X_o_cutoff_*(IMS_bfCut_ + IMS_dfCut_ * IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_efCut_ = - eterm * (tmp);
        if (fabs(IMS_efCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("IdealMolalSoln::calcCutoffParams_",
                           "failed to converge on the f polynomial");
    }
    converged = false;
    double f_0 = IMS_afCut_ + IMS_efCut_;
    double f_prime_0 = 1.0 - IMS_afCut_ / IMS_cCut_ + IMS_bfCut_;
    IMS_egCut_ = 0.0;
    for (int its = 0; its < 100 && !converged; its++) {
        double oldV = IMS_egCut_;
        double lng_0 = -log(IMS_gamma_o_min_) - f_prime_0 / f_0;
        IMS_agCut_ = exp(lng_0) - IMS_egCut_;
        IMS_bgCut_ = IMS_agCut_ / IMS_cCut_ + IMS_slopegCut_ - 1.0;
        IMS_dgCut_ = ((- IMS_agCut_/IMS_cCut_ + IMS_bgCut_ - IMS_bgCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_agCut_ + IMS_X_o_cutoff_*(IMS_bgCut_ + IMS_dgCut_ *IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_egCut_ = - eterm * (tmp);
        if (fabs(IMS_egCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("IdealMolalSoln::calcCutoffParams_",
                           "failed to converge on the g polynomial");
    }
}

}
