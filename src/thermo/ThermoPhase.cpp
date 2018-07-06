/**
 *  @file ThermoPhase.cpp
 * Definition file for class ThermoPhase, the base class for phases with
 * thermodynamic properties
 * (see class \link Cantera::ThermoPhase ThermoPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/equil/ChemEquil.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/base/ctml.h"

#include <iomanip>
#include <fstream>

using namespace std;

namespace Cantera
{

ThermoPhase::ThermoPhase() :
    m_speciesData(0),
    m_phi(0.0),
    m_hasElementPotentials(false),
    m_chargeNeutralityNecessary(false),
    m_ssConvention(cSS_CONVENTION_TEMPERATURE),
    m_tlast(0.0)
{
}

ThermoPhase::~ThermoPhase()
{
    for (size_t k = 0; k < m_speciesData.size(); k++) {
        delete m_speciesData[k];
    }
}

void ThermoPhase::resetHf298(size_t k) {
    if (k != npos) {
        m_spthermo.resetHf298(k);
    } else {
        for (size_t k = 0; k < nSpecies(); k++) {
            m_spthermo.resetHf298(k);
        }
    }
    invalidateCache();
}

int ThermoPhase::activityConvention() const
{
    return cAC_CONVENTION_MOLAR;
}

int ThermoPhase::standardStateConvention() const
{
    return m_ssConvention;
}

doublereal ThermoPhase::logStandardConc(size_t k) const
{
    return log(standardConcentration(k));
}

void ThermoPhase::getActivities(doublereal* a) const
{
    getActivityConcentrations(a);
    for (size_t k = 0; k < nSpecies(); k++) {
        a[k] /= standardConcentration(k);
    }
}

void ThermoPhase::getLnActivityCoefficients(doublereal* lnac) const
{
    getActivityCoefficients(lnac);
    for (size_t k = 0; k < m_kk; k++) {
        lnac[k] = std::log(lnac[k]);
    }
}

void ThermoPhase::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}

void ThermoPhase::setState_TPX(doublereal t, doublereal p, const doublereal* x)
{
    setMoleFractions(x);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPX(doublereal t, doublereal p, const compositionMap& x)
{
    setMoleFractionsByName(x);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPX(doublereal t, doublereal p, const std::string& x)
{
    setMoleFractionsByName(x);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPY(doublereal t, doublereal p, const doublereal* y)
{
    setMassFractions(y);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPY(doublereal t, doublereal p, const compositionMap& y)
{
    setMassFractionsByName(y);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPY(doublereal t, doublereal p, const std::string& y)
{
    setMassFractionsByName(y);
    setState_TP(t,p);
}

void ThermoPhase::setState_TP(doublereal t, doublereal p)
{
    setTemperature(t);
    setPressure(p);
}

void ThermoPhase::setState_RPX(doublereal rho, doublereal p, const doublereal* x)
{
    setMoleFractions(x);
    setState_RP(rho, p);
}

void ThermoPhase::setState_RPX(doublereal rho, doublereal p, const compositionMap& x)
{
    setMoleFractionsByName(x);
    setState_RP(rho,p);
}

void ThermoPhase::setState_RPX(doublereal rho, doublereal p, const std::string& x)
{
    setMoleFractionsByName(x);
    setState_RP(rho,p);
}

void ThermoPhase::setState_RPY(doublereal rho, doublereal p, const doublereal* y)
{
    setMassFractions(y);
    setState_RP(rho,p);
}

void ThermoPhase::setState_RPY(doublereal rho, doublereal p, const compositionMap& y)
{
    setMassFractionsByName(y);
    setState_RP(rho,p);
}

void ThermoPhase::setState_RPY(doublereal rho, doublereal p, const std::string& y)
{
    setMassFractionsByName(y);
    setState_RP(rho,p);
}

void ThermoPhase::setState_PX(doublereal p, doublereal* x)
{
    setMoleFractions(x);
    setPressure(p);
}

void ThermoPhase::setState_PY(doublereal p, doublereal* y)
{
    setMassFractions(y);
    setPressure(p);
}

void ThermoPhase::setState_HP(double Htarget, double p, double rtol)
{
    setState_HPorUV(Htarget, p, rtol, false);
}

void ThermoPhase::setState_UV(double u, double v, double rtol)
{
    setState_HPorUV(u, v, rtol, true);
}

void ThermoPhase::setState_conditional_TP(doublereal t, doublereal p, bool set_p)
{
    setTemperature(t);
    if (set_p) {
        setPressure(p);
    }
}

void ThermoPhase::setState_HPorUV(double Htarget, double p,
                                  double rtol, bool doUV)
{
    doublereal dt;
    doublereal v = 0.0;

    // Assign the specific volume or pressure and make sure it's positive
    if (doUV) {
        doublereal v = p;
        if (v < 1.0E-300) {
            throw CanteraError("setState_HPorUV (UV)",
                               "Input specific volume is too small or negative. v = {}", v);
        }
        setDensity(1.0/v);
    } else {
        if (p < 1.0E-300) {
            throw CanteraError("setState_HPorUV (HP)",
                               "Input pressure is too small or negative. p = {}", p);
        }
        setPressure(p);
    }
    double Tmax = maxTemp() + 0.1;
    double Tmin = minTemp() - 0.1;

    // Make sure we are within the temperature bounds at the start
    // of the iteration
    double Tnew = temperature();
    double Tinit = Tnew;
    if (Tnew > Tmax) {
        Tnew = Tmax - 1.0;
    } else if (Tnew < Tmin) {
        Tnew = Tmin + 1.0;
    }
    if (Tnew != Tinit) {
        setState_conditional_TP(Tnew, p, !doUV);
    }

    double Hnew = (doUV) ? intEnergy_mass() : enthalpy_mass();
    double Cpnew = (doUV) ? cv_mass() : cp_mass();
    double Htop = Hnew;
    double Ttop = Tnew;
    double Hbot = Hnew;
    double Tbot = Tnew;

    bool ignoreBounds = false;
    // Unstable phases are those for which cp < 0.0. These are possible for
    // cases where we have passed the spinodal curve.
    bool unstablePhase = false;
    // Counter indicating the last temperature point where the
    // phase was unstable
    double Tunstable = -1.0;
    bool unstablePhaseNew = false;

    // Newton iteration
    for (int n = 0; n < 500; n++) {
        double Told = Tnew;
        double Hold = Hnew;
        double cpd = Cpnew;
        if (cpd < 0.0) {
            unstablePhase = true;
            Tunstable = Tnew;
        }
        // limit step size to 100 K
        dt = clip((Htarget - Hold)/cpd, -100.0, 100.0);

        // Calculate the new T
        Tnew = Told + dt;

        // Limit the step size so that we are convergent This is the step that
        // makes it different from a Newton's algorithm
        if ((dt > 0.0 && unstablePhase) || (dt <= 0.0 && !unstablePhase)) {
            if (Hbot < Htarget && Tnew < (0.75 * Tbot + 0.25 * Told)) {
                dt = 0.75 * (Tbot - Told);
                Tnew = Told + dt;
            }
        } else if (Htop > Htarget && Tnew > (0.75 * Ttop + 0.25 * Told)) {
            dt = 0.75 * (Ttop - Told);
            Tnew = Told + dt;
        }

        // Check Max and Min values
        if (Tnew > Tmax && !ignoreBounds) {
            setState_conditional_TP(Tmax, p, !doUV);
            double Hmax = (doUV) ? intEnergy_mass() : enthalpy_mass();
            if (Hmax >= Htarget) {
                if (Htop < Htarget) {
                    Ttop = Tmax;
                    Htop = Hmax;
                }
            } else {
                Tnew = Tmax + 1.0;
                ignoreBounds = true;
            }
        }
        if (Tnew < Tmin && !ignoreBounds) {
            setState_conditional_TP(Tmin, p, !doUV);
            double Hmin = (doUV) ? intEnergy_mass() : enthalpy_mass();
            if (Hmin <= Htarget) {
                if (Hbot > Htarget) {
                    Tbot = Tmin;
                    Hbot = Hmin;
                }
            } else {
                Tnew = Tmin - 1.0;
                ignoreBounds = true;
            }
        }

        // Try to keep phase within its region of stability
        // -> Could do a lot better if I calculate the
        //    spinodal value of H.
        for (int its = 0; its < 10; its++) {
            Tnew = Told + dt;
            if (Tnew < Told / 3.0) {
                Tnew = Told / 3.0;
                dt = -2.0 * Told / 3.0;
            }
            setState_conditional_TP(Tnew, p, !doUV);
            if (doUV) {
                Hnew = intEnergy_mass();
                Cpnew = cv_mass();
            } else {
                Hnew = enthalpy_mass();
                Cpnew = cp_mass();
            }
            if (Cpnew < 0.0) {
                unstablePhaseNew = true;
                Tunstable = Tnew;
            } else {
                unstablePhaseNew = false;
                break;
            }
            if (unstablePhase == false && unstablePhaseNew == true) {
                dt *= 0.25;
            }
        }

        if (Hnew == Htarget) {
            return;
        } else if (Hnew > Htarget && (Htop < Htarget || Hnew < Htop)) {
            Htop = Hnew;
            Ttop = Tnew;
        } else if (Hnew < Htarget && (Hbot > Htarget || Hnew > Hbot)) {
            Hbot = Hnew;
            Tbot = Tnew;
        }
        // Convergence in H
        double Herr = Htarget - Hnew;
        double acpd = std::max(fabs(cpd), 1.0E-5);
        double denom = std::max(fabs(Htarget), acpd * Tnew);
        double HConvErr = fabs((Herr)/denom);
        if (HConvErr < rtol || fabs(dt/Tnew) < rtol) {
            return;
        }
    }
    // We are here when there hasn't been convergence

    // Formulate a detailed error message, since questions seem to arise often
    // about the lack of convergence.
    string ErrString =  "No convergence in 500 iterations\n";
    if (doUV) {
        ErrString += fmt::format(
            "\tTarget Internal Energy  = {}\n"
            "\tCurrent Specific Volume = {}\n"
            "\tStarting Temperature    = {}\n"
            "\tCurrent Temperature     = {}\n"
            "\tCurrent Internal Energy = {}\n"
            "\tCurrent Delta T         = {}\n",
            Htarget, v, Tinit, Tnew, Hnew, dt);
    } else {
        ErrString += fmt::format(
            "\tTarget Enthalpy         = {}\n"
            "\tCurrent Pressure        = {}\n"
            "\tStarting Temperature    = {}\n"
            "\tCurrent Temperature     = {}\n"
            "\tCurrent Enthalpy        = {}\n"
            "\tCurrent Delta T         = {}\n",
            Htarget, p, Tinit, Tnew, Hnew, dt);
    }
    if (unstablePhase) {
        ErrString += fmt::format(
            "\t  - The phase became unstable (Cp < 0) T_unstable_last = {}\n",
            Tunstable);
    }
    if (doUV) {
        throw CanteraError("setState_HPorUV (UV)", ErrString);
    } else {
        throw CanteraError("setState_HPorUV (HP)", ErrString);
    }
}

void ThermoPhase::setState_SP(double Starget, double p, double rtol)
{
    setState_SPorSV(Starget, p, rtol, false);
}

void ThermoPhase::setState_SV(double Starget, double v, double rtol)
{
    setState_SPorSV(Starget, v, rtol, true);
}

void ThermoPhase::setState_SPorSV(double Starget, double p,
                                  double rtol, bool doSV)
{
    doublereal v = 0.0;
    doublereal dt;
    if (doSV) {
        v = p;
        if (v < 1.0E-300) {
            throw CanteraError("setState_SPorSV (SV)",
                "Input specific volume is too small or negative. v = {}", v);
        }
        setDensity(1.0/v);
    } else {
        if (p < 1.0E-300) {
            throw CanteraError("setState_SPorSV (SP)",
                "Input pressure is too small or negative. p = {}", p);
        }
        setPressure(p);
    }
    double Tmax = maxTemp() + 0.1;
    double Tmin = minTemp() - 0.1;

    // Make sure we are within the temperature bounds at the start
    // of the iteration
    double Tnew = temperature();
    double Tinit = Tnew;
    if (Tnew > Tmax) {
        Tnew = Tmax - 1.0;
    } else if (Tnew < Tmin) {
        Tnew = Tmin + 1.0;
    }
    if (Tnew != Tinit) {
        setState_conditional_TP(Tnew, p, !doSV);
    }

    double Snew = entropy_mass();
    double Cpnew = (doSV) ? cv_mass() : cp_mass();
    double Stop = Snew;
    double Ttop = Tnew;
    double Sbot = Snew;
    double Tbot = Tnew;

    bool ignoreBounds = false;
    // Unstable phases are those for which Cp < 0.0. These are possible for
    // cases where we have passed the spinodal curve.
    bool unstablePhase = false;
    double Tunstable = -1.0;
    bool unstablePhaseNew = false;

    // Newton iteration
    for (int n = 0; n < 500; n++) {
        double Told = Tnew;
        double Sold = Snew;
        double cpd = Cpnew;
        if (cpd < 0.0) {
            unstablePhase = true;
            Tunstable = Tnew;
        }
        // limit step size to 100 K
        dt = clip((Starget - Sold)*Told/cpd, -100.0, 100.0);
        Tnew = Told + dt;

        // Limit the step size so that we are convergent
        if ((dt > 0.0 && unstablePhase) || (dt <= 0.0 && !unstablePhase)) {
            if (Sbot < Starget && Tnew < Tbot) {
                dt = 0.75 * (Tbot - Told);
                Tnew = Told + dt;
            }
        } else if (Stop > Starget && Tnew > Ttop) {
            dt = 0.75 * (Ttop - Told);
            Tnew = Told + dt;
        }

        // Check Max and Min values
        if (Tnew > Tmax && !ignoreBounds) {
            setState_conditional_TP(Tmax, p, !doSV);
            double Smax = entropy_mass();
            if (Smax >= Starget) {
                if (Stop < Starget) {
                    Ttop = Tmax;
                    Stop = Smax;
                }
            } else {
                Tnew = Tmax + 1.0;
                ignoreBounds = true;
            }
        } else if (Tnew < Tmin && !ignoreBounds) {
            setState_conditional_TP(Tmin, p, !doSV);
            double Smin = entropy_mass();
            if (Smin <= Starget) {
                if (Sbot > Starget) {
                    Tbot = Tmin;
                    Sbot = Smin;
                }
            } else {
                Tnew = Tmin - 1.0;
                ignoreBounds = true;
            }
        }

        // Try to keep phase within its region of stability
        // -> Could do a lot better if I calculate the
        //    spinodal value of H.
        for (int its = 0; its < 10; its++) {
            Tnew = Told + dt;
            setState_conditional_TP(Tnew, p, !doSV);
            Cpnew = (doSV) ? cv_mass() : cp_mass();
            Snew = entropy_mass();
            if (Cpnew < 0.0) {
                unstablePhaseNew = true;
                Tunstable = Tnew;
            } else {
                unstablePhaseNew = false;
                break;
            }
            if (unstablePhase == false && unstablePhaseNew == true) {
                dt *= 0.25;
            }
        }

        if (Snew == Starget) {
            return;
        } else if (Snew > Starget && (Stop < Starget || Snew < Stop)) {
            Stop = Snew;
            Ttop = Tnew;
        } else if (Snew < Starget && (Sbot > Starget || Snew > Sbot)) {
            Sbot = Snew;
            Tbot = Tnew;
        }
        // Convergence in S
        double Serr = Starget - Snew;
        double acpd = std::max(fabs(cpd), 1.0E-5);
        double denom = std::max(fabs(Starget), acpd * Tnew);
        double SConvErr = fabs((Serr * Tnew)/denom);
        if (SConvErr < rtol || fabs(dt/Tnew) < rtol) {
            return;
        }
    }
    // We are here when there hasn't been convergence

    // Formulate a detailed error message, since questions seem to arise often
    // about the lack of convergence.
    string ErrString =  "No convergence in 500 iterations\n";
    if (doSV) {
        ErrString += fmt::format(
            "\tTarget Entropy          = {}\n"
            "\tCurrent Specific Volume = {}\n"
            "\tStarting Temperature    = {}\n"
            "\tCurrent Temperature     = {}\n"
            "\tCurrent Entropy         = {}\n"
            "\tCurrent Delta T         = {}\n",
            Starget, v, Tinit, Tnew, Snew, dt);
    } else {
        ErrString += fmt::format(
            "\tTarget Entropy          = {}\n"
            "\tCurrent Pressure        = {}\n"
            "\tStarting Temperature    = {}\n"
            "\tCurrent Temperature     = {}\n"
            "\tCurrent Entropy         = {}\n"
            "\tCurrent Delta T         = {}\n",
            Starget, p, Tinit, Tnew, Snew, dt);
    }
    if (unstablePhase) {
        ErrString += fmt::format("\t  - The phase became unstable (Cp < 0) T_unstable_last = {}\n",
                     Tunstable);
    }
    if (doSV) {
        throw CanteraError("setState_SPorSV (SV)", ErrString);
    } else {
        throw CanteraError("setState_SPorSV (SP)", ErrString);
    }
}

MultiSpeciesThermo& ThermoPhase::speciesThermo(int k)
{
    return m_spthermo;
}

void ThermoPhase::initThermoFile(const std::string& inputFile,
                                 const std::string& id)
{
    XML_Node* fxml = get_XML_File(inputFile);
    XML_Node* fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
        throw CanteraError("ThermoPhase::initThermoFile",
                           "ERROR: Can not find phase named {} in file"
                           " named {}", id, inputFile);
    }
    importPhase(*fxml_phase, this);
}

void ThermoPhase::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    if (phaseNode.hasChild("state")) {
        setStateFromXML(phaseNode.child("state"));
    }
}

void ThermoPhase::initThermo()
{
    // Check to see that all of the species thermo objects have been initialized
    if (!m_spthermo.ready(m_kk)) {
        throw CanteraError("ThermoPhase::initThermo()",
                           "Missing species thermo data");
    }
}

bool ThermoPhase::addSpecies(shared_ptr<Species> spec)
{
    if (!spec->thermo) {
        throw CanteraError("ThermoPhase::addSpecies",
            "Species {} has no thermo data", spec->name);
    }
    bool added = Phase::addSpecies(spec);
    if (added) {
        spec->thermo->validate(spec->name);
        m_spthermo.install_STIT(m_kk-1, spec->thermo);
    }
    return added;
}

void ThermoPhase::modifySpecies(size_t k, shared_ptr<Species> spec)
{
    if (!spec->thermo) {
        throw CanteraError("ThermoPhase::modifySpecies",
            "Species {} has no thermo data", spec->name);
    }
    Phase::modifySpecies(k, spec);
    if (speciesName(k) != spec->name) {
        throw CanteraError("ThermoPhase::modifySpecies",
            "New species '{}' does not match existing species '{}' at index {}",
                           spec->name, speciesName(k), k);
    }
    spec->thermo->validate(spec->name);
    m_spthermo.modifySpecies(k, spec->thermo);
}

void ThermoPhase::saveSpeciesData(const size_t k, const XML_Node* const data)
{
    if (m_speciesData.size() < (k + 1)) {
        m_speciesData.resize(k+1, 0);
    }
    m_speciesData[k] = new XML_Node(*data);
}

const std::vector<const XML_Node*> & ThermoPhase::speciesData() const
{
    if (m_speciesData.size() != m_kk) {
        throw CanteraError("ThermoPhase::speciesData",
                           "m_speciesData is the wrong size");
    }
    return m_speciesData;
}

void ThermoPhase::setStateFromXML(const XML_Node& state)
{
    string comp = getChildValue(state,"moleFractions");
    if (comp != "") {
        setMoleFractionsByName(comp);
    } else {
        comp = getChildValue(state,"massFractions");
        if (comp != "") {
            setMassFractionsByName(comp);
        }
    }
    if (state.hasChild("temperature")) {
        double t = getFloat(state, "temperature", "temperature");
        setTemperature(t);
    }
    if (state.hasChild("pressure")) {
        double p = getFloat(state, "pressure", "pressure");
        setPressure(p);
    }
    if (state.hasChild("density")) {
        double rho = getFloat(state, "density", "density");
        setDensity(rho);
    }
}

void ThermoPhase::invalidateCache() {
    Phase::invalidateCache();
    m_tlast += 0.1234;
}

void ThermoPhase::equilibrate(const std::string& XY, const std::string& solver,
                              double rtol, int max_steps, int max_iter,
                              int estimate_equil, int log_level)
{
    if (solver == "auto" || solver == "element_potential") {
        vector_fp initial_state;
        saveState(initial_state);
        debuglog("Trying ChemEquil solver\n", log_level);
        try {
            ChemEquil E;
            E.options.maxIterations = max_steps;
            E.options.relTolerance = rtol;
            int ret = E.equilibrate(*this, XY.c_str(), log_level-1);
            if (ret < 0) {
                throw CanteraError("ThermoPhase::equilibrate",
                    "ChemEquil solver failed. Return code: {}", ret);
            }
            m_lambdaRRT.resize(nElements());
            for (size_t m = 0; m < nElements(); m++) {
                m_lambdaRRT[m] = E.elementPotentials()[m] / RT();
            }
            m_hasElementPotentials = true;
            debuglog("ChemEquil solver succeeded\n", log_level);
            return;
        } catch (std::exception& err) {
            debuglog("ChemEquil solver failed.\n", log_level);
            debuglog(err.what(), log_level);
            restoreState(initial_state);
            if (solver == "auto") {
            } else {
                throw;
            }
        }
    }

    if (solver == "auto" || solver == "vcs" || solver == "gibbs") {
        MultiPhase M;
        M.addPhase(this, 1.0);
        M.init();
        M.equilibrate(XY, solver, rtol, max_steps, max_iter,
                      estimate_equil, log_level);
        return;
    }

    if (solver != "auto") {
        throw CanteraError("ThermoPhase::equilibrate",
                           "Invalid solver specified: '{}'", solver);
    }
}

void ThermoPhase::setElementPotentials(const vector_fp& lambda)
{
    warn_deprecated("ThermoPhase::setElementPotentials",
        "To be removed after Cantera 2.4");
    size_t mm = nElements();
    if (lambda.size() < mm) {
        throw CanteraError("setElementPotentials", "lambda too small");
    }
    if (!m_hasElementPotentials) {
        m_lambdaRRT.resize(mm);
    }
    scale(lambda.begin(), lambda.end(), m_lambdaRRT.begin(), 1.0/RT());
    m_hasElementPotentials = true;
}

bool ThermoPhase::getElementPotentials(doublereal* lambda) const
{
    warn_deprecated("ThermoPhase::getElementPotentials",
        "To be removed after Cantera 2.4");
    if (m_hasElementPotentials) {
        scale(m_lambdaRRT.begin(), m_lambdaRRT.end(), lambda, RT());
    }
    return m_hasElementPotentials;
}

void ThermoPhase::getdlnActCoeffdlnN(const size_t ld, doublereal* const dlnActCoeffdlnN)
{
    for (size_t m = 0; m < m_kk; m++) {
        for (size_t k = 0; k < m_kk; k++) {
            dlnActCoeffdlnN[ld * k + m] = 0.0;
        }
    }
    return;
}

void ThermoPhase::getdlnActCoeffdlnN_numderiv(const size_t ld, doublereal* const dlnActCoeffdlnN)
{
    double deltaMoles_j = 0.0;
    double pres = pressure();

    // Evaluate the current base activity coefficients if necessary
    vector_fp ActCoeff_Base(m_kk);
    getActivityCoefficients(ActCoeff_Base.data());
    vector_fp Xmol_Base(m_kk);
    getMoleFractions(Xmol_Base.data());

    // Make copies of ActCoeff and Xmol_ for use in taking differences
    vector_fp ActCoeff(m_kk);
    vector_fp Xmol(m_kk);
    double v_totalMoles = 1.0;
    double TMoles_base = v_totalMoles;

    // Loop over the columns species to be deltad
    for (size_t j = 0; j < m_kk; j++) {
        // Calculate a value for the delta moles of species j
        // -> Note Xmol_[] and Tmoles are always positive or zero quantities.
        // -> experience has shown that you always need to make the deltas
        //    greater than needed to change the other mole fractions in order
        //    to capture some effects.
        double moles_j_base = v_totalMoles * Xmol_Base[j];
        deltaMoles_j = 1.0E-7 * moles_j_base + v_totalMoles * 1.0E-13 + 1.0E-150;

        // Now, update the total moles in the phase and all of the mole
        // fractions based on this.
        v_totalMoles = TMoles_base + deltaMoles_j;
        for (size_t k = 0; k < m_kk; k++) {
            Xmol[k] = Xmol_Base[k] * TMoles_base / v_totalMoles;
        }
        Xmol[j] = (moles_j_base + deltaMoles_j) / v_totalMoles;

        // Go get new values for the activity coefficients.
        // -> Note this calls setState_PX();
        setState_PX(pres, Xmol.data());
        getActivityCoefficients(ActCoeff.data());

        // Calculate the column of the matrix
        double* const lnActCoeffCol = dlnActCoeffdlnN + ld * j;
        for (size_t k = 0; k < m_kk; k++) {
            lnActCoeffCol[k] = (2*moles_j_base + deltaMoles_j) *(ActCoeff[k] - ActCoeff_Base[k]) /
                               ((ActCoeff[k] + ActCoeff_Base[k]) * deltaMoles_j);
        }
        // Revert to the base case Xmol_, v_totalMoles
        v_totalMoles = TMoles_base;
        Xmol = Xmol_Base;
    }

    setState_PX(pres, Xmol_Base.data());
}

std::string ThermoPhase::report(bool show_thermo, doublereal threshold) const
{
    fmt::memory_buffer b;
    try {
        if (name() != "") {
            format_to(b, "\n  {}:\n", name());
        }
        format_to(b, "\n");
        format_to(b, "       temperature    {:12.6g}  K\n", temperature());
        format_to(b, "          pressure    {:12.6g}  Pa\n", pressure());
        format_to(b, "           density    {:12.6g}  kg/m^3\n", density());
        format_to(b, "  mean mol. weight    {:12.6g}  amu\n", meanMolecularWeight());

        doublereal phi = electricPotential();
        if (phi != 0.0) {
            format_to(b, "         potential    {:12.6g}  V\n", phi);
        }
        if (show_thermo) {
            format_to(b, "\n");
            format_to(b, "                          1 kg            1 kmol\n");
            format_to(b, "                       -----------      ------------\n");
            format_to(b, "          enthalpy    {:12.5g}     {:12.4g}     J\n",
                    enthalpy_mass(), enthalpy_mole());
            format_to(b, "   internal energy    {:12.5g}     {:12.4g}     J\n",
                    intEnergy_mass(), intEnergy_mole());
            format_to(b, "           entropy    {:12.5g}     {:12.4g}     J/K\n",
                    entropy_mass(), entropy_mole());
            format_to(b, "    Gibbs function    {:12.5g}     {:12.4g}     J\n",
                    gibbs_mass(), gibbs_mole());
            format_to(b, " heat capacity c_p    {:12.5g}     {:12.4g}     J/K\n",
                    cp_mass(), cp_mole());
            try {
                format_to(b, " heat capacity c_v    {:12.5g}     {:12.4g}     J/K\n",
                        cv_mass(), cv_mole());
            } catch (NotImplementedError&) {
                format_to(b, " heat capacity c_v    <not implemented>       \n");
            }
        }

        vector_fp x(m_kk);
        vector_fp y(m_kk);
        vector_fp mu(m_kk);
        getMoleFractions(&x[0]);
        getMassFractions(&y[0]);
        getChemPotentials(&mu[0]);
        int nMinor = 0;
        doublereal xMinor = 0.0;
        doublereal yMinor = 0.0;
        format_to(b, "\n");
        if (show_thermo) {
            format_to(b, "                           X     "
                    "            Y          Chem. Pot. / RT\n");
            format_to(b, "                     -------------     "
                    "------------     ------------\n");
            for (size_t k = 0; k < m_kk; k++) {
                if (abs(x[k]) >= threshold) {
                    if (abs(x[k]) > SmallNumber) {
                        format_to(b, "{:>18s}   {:12.6g}     {:12.6g}     {:12.6g}\n",
                                speciesName(k), x[k], y[k], mu[k]/RT());
                    } else {
                        format_to(b, "{:>18s}   {:12.6g}     {:12.6g}\n",
                                speciesName(k), x[k], y[k]);
                    }
                } else {
                    nMinor++;
                    xMinor += x[k];
                    yMinor += y[k];
                }
            }
        } else {
            format_to(b, "                           X                 Y\n");
            format_to(b, "                     -------------     ------------\n");
            for (size_t k = 0; k < m_kk; k++) {
                if (abs(x[k]) >= threshold) {
                    format_to(b, "{:>18s}   {:12.6g}     {:12.6g}\n",
                            speciesName(k), x[k], y[k]);
                } else {
                    nMinor++;
                    xMinor += x[k];
                    yMinor += y[k];
                }
            }
        }
        if (nMinor) {
            format_to(b, "     [{:+5d} minor]   {:12.6g}     {:12.6g}\n",
                    nMinor, xMinor, yMinor);
        }
    } catch (CanteraError& err) {
        return to_string(b) + err.what();
    }
    return to_string(b);
}

void ThermoPhase::reportCSV(std::ofstream& csvFile) const
{
    int tabS = 15;
    int tabM = 30;
    csvFile.precision(8);
    vector_fp X(nSpecies());
    getMoleFractions(&X[0]);
    std::vector<std::string> pNames;
    std::vector<vector_fp> data;
    getCsvReportData(pNames, data);

    csvFile << setw(tabS) << "Species,";
    for (size_t i = 0; i < pNames.size(); i++) {
        csvFile << setw(tabM) << pNames[i] << ",";
    }
    csvFile << endl;
    for (size_t k = 0; k < nSpecies(); k++) {
        csvFile << setw(tabS) << speciesName(k) + ",";
        if (X[k] > SmallNumber) {
            for (size_t i = 0; i < pNames.size(); i++) {
                csvFile << setw(tabM) << data[i][k] << ",";
            }
            csvFile << endl;
        } else {
            for (size_t i = 0; i < pNames.size(); i++) {
                csvFile << setw(tabM) << 0 << ",";
            }
            csvFile << endl;
        }
    }
}

void ThermoPhase::getCsvReportData(std::vector<std::string>& names,
                                   std::vector<vector_fp>& data) const
{
    names.clear();
    data.assign(10, vector_fp(nSpecies()));

    names.push_back("X");
    getMoleFractions(&data[0][0]);

    names.push_back("Y");
    getMassFractions(&data[1][0]);

    names.push_back("Chem. Pot (J/kmol)");
    getChemPotentials(&data[2][0]);

    names.push_back("Activity");
    getActivities(&data[3][0]);

    names.push_back("Act. Coeff.");
    getActivityCoefficients(&data[4][0]);

    names.push_back("Part. Mol Enthalpy (J/kmol)");
    getPartialMolarEnthalpies(&data[5][0]);

    names.push_back("Part. Mol. Entropy (J/K/kmol)");
    getPartialMolarEntropies(&data[6][0]);

    names.push_back("Part. Mol. Energy (J/kmol)");
    getPartialMolarIntEnergies(&data[7][0]);

    names.push_back("Part. Mol. Cp (J/K/kmol");
    getPartialMolarCp(&data[8][0]);

    names.push_back("Part. Mol. Cv (J/K/kmol)");
    getPartialMolarVolumes(&data[9][0]);
}

}
