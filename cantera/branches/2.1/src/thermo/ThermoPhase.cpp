/**
 *  @file ThermoPhase.cpp
 * Definition file for class ThermoPhase, the base class for phases with
 * thermodynamic properties
 * (see class \link Cantera::ThermoPhase ThermoPhase\endlink).
 */

//  Copyright 2002 California Institute of Technology

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoFactory.h"

#include <iomanip>
#include <fstream>

using namespace std;
using namespace ctml;

namespace Cantera
{

ThermoPhase::ThermoPhase() :
    Phase(),
    m_spthermo(0), m_speciesData(0),
    m_phi(0.0),
    m_hasElementPotentials(false),
    m_chargeNeutralityNecessary(false),
    m_ssConvention(cSS_CONVENTION_TEMPERATURE)
{
}

ThermoPhase::~ThermoPhase()
{
    for (size_t k = 0; k < m_speciesData.size(); k++) {
        delete m_speciesData[k];
    }
    delete m_spthermo;
}

ThermoPhase::ThermoPhase(const ThermoPhase& right)  :
    Phase(),
    m_spthermo(0),
    m_speciesData(0),
    m_phi(0.0),
    m_hasElementPotentials(false),
    m_chargeNeutralityNecessary(false),
    m_ssConvention(cSS_CONVENTION_TEMPERATURE)
{
    /*
     * Call the assignment operator
     */
    *this = operator=(right);
}

ThermoPhase& ThermoPhase::
operator=(const ThermoPhase& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    /*
     * We need to destruct first
     */
    for (size_t k = 0; k < m_speciesData.size(); k++) {
        delete m_speciesData[k];
    }
    delete m_spthermo;

    /*
     * Call the base class assignment operator
     */
    (void)Phase::operator=(right);

    /*
     * Pointer to the species thermodynamic property manager
     * We own this, so we need to do a deep copy
     */
    m_spthermo = (right.m_spthermo)->duplMyselfAsSpeciesThermo();

    /*
     * Do a deep copy of species Data, because we own this
     */
    m_speciesData.resize(m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        m_speciesData[k] = new XML_Node(*(right.m_speciesData[k]));
    }

    m_phi = right.m_phi;
    m_lambdaRRT = right.m_lambdaRRT;
    m_hasElementPotentials = right.m_hasElementPotentials;
    m_chargeNeutralityNecessary = right.m_chargeNeutralityNecessary;
    m_ssConvention = right.m_ssConvention;
    return *this;
}

ThermoPhase* ThermoPhase::duplMyselfAsThermoPhase() const
{
    return new ThermoPhase(*this);
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

void ThermoPhase::setState_TPX(doublereal t, doublereal p, const doublereal* x)
{
    setMoleFractions(x);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPX(doublereal t, doublereal p, compositionMap& x)
{
    setMoleFractionsByName(x);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPX(doublereal t, doublereal p, const std::string& x)
{
    compositionMap xx = parseCompString(x, speciesNames());
    setMoleFractionsByName(xx);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPY(doublereal t, doublereal p,
                               const doublereal* y)
{
    setMassFractions(y);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPY(doublereal t, doublereal p,
                               compositionMap& y)
{
    setMassFractionsByName(y);
    setState_TP(t,p);
}

void ThermoPhase::setState_TPY(doublereal t, doublereal p,
                               const std::string& y)
{
    compositionMap yy = parseCompString(y, speciesNames());
    setMassFractionsByName(yy);
    setState_TP(t,p);
}

void ThermoPhase::setState_TP(doublereal t, doublereal p)
{
    setTemperature(t);
    setPressure(p);
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

void ThermoPhase::setState_HP(doublereal Htarget, doublereal p,
                              doublereal dTtol)
{
    setState_HPorUV(Htarget, p, dTtol, false);
}

void ThermoPhase::setState_UV(doublereal u, doublereal v,
                              doublereal dTtol)
{
    setState_HPorUV(u, v, dTtol, true);
}

void ThermoPhase::setState_conditional_TP(doublereal t, doublereal p, bool set_p)
{
    setTemperature(t);
    if (set_p) {
        setPressure(p);
    }
}

void ThermoPhase::setState_HPorUV(doublereal Htarget, doublereal p,
                                  doublereal dTtol, bool doUV)
{
    doublereal dt;
    doublereal Hmax = 0.0, Hmin = 0.0;
    doublereal v = 0.0;

    // Assign the specific volume or pressure and make sure it's positive
    if (doUV) {
        v = p;
        if (v < 1.0E-300) {
            throw CanteraError("setState_HPorUV (UV)",
                               "Input specific volume is too small or negative. v = " + fp2str(v));
        }
        setDensity(1.0/v);
    } else {
        if (p < 1.0E-300) {
            throw CanteraError("setState_HPorUV (HP)",
                               "Input pressure is too small or negative. p = " + fp2str(p));
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
    double Told = Tnew;
    double Hold = Hnew;

    bool ignoreBounds = false;
    // Unstable phases are those for which
    // cp < 0.0. These are possible for cases where
    // we have passed the spinodal curve.
    bool unstablePhase = false;
    // Counter indicating the last temperature point where the
    // phase was unstable
    double Tunstable = -1.0;
    bool unstablePhaseNew = false;

    // Newton iteration
    for (int n = 0; n < 500; n++) {
        Told = Tnew;
        Hold = Hnew;
        double cpd = Cpnew;
        if (cpd < 0.0) {
            unstablePhase = true;
            Tunstable = Tnew;
        }
        // limit step size to 100 K
        dt = clip((Htarget - Hold)/cpd, -100.0, 100.0);

        // Calculate the new T
        Tnew = Told + dt;

        // Limit the step size so that we are convergent
        // This is the step that makes it different from a
        // Newton's algorithm
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
            Hmax = (doUV) ? intEnergy_mass() : enthalpy_mass();
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
            Hmin = (doUV) ? intEnergy_mass() : enthalpy_mass();
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
        double denom = std::max(fabs(Htarget), acpd * dTtol);
        double HConvErr = fabs((Herr)/denom);
        if (HConvErr < 0.00001 *dTtol || fabs(dt) < dTtol) {
            return;
        }
    }
    // We are here when there hasn't been convergence
    /*
     * Formulate a detailed error message, since questions seem to
     * arise often about the lack of convergence.
     */
    string ErrString =  "No convergence in 500 iterations\n";
    if (doUV) {
        ErrString += "\tTarget Internal Energy  = " + fp2str(Htarget) + "\n";
        ErrString += "\tCurrent Specific Volume = " + fp2str(v) + "\n";
        ErrString += "\tStarting Temperature    = " + fp2str(Tinit) + "\n";
        ErrString += "\tCurrent Temperature     = " + fp2str(Tnew) + "\n";
        ErrString += "\tCurrent Internal Energy = " + fp2str(Hnew) + "\n";
        ErrString += "\tCurrent Delta T         = " + fp2str(dt) + "\n";
    } else {
        ErrString += "\tTarget Enthalpy         = " + fp2str(Htarget) + "\n";
        ErrString += "\tCurrent Pressure        = " + fp2str(p) + "\n";
        ErrString += "\tStarting Temperature    = " + fp2str(Tinit) + "\n";
        ErrString += "\tCurrent Temperature     = " + fp2str(Tnew) + "\n";
        ErrString += "\tCurrent Enthalpy        = " + fp2str(Hnew) + "\n";
        ErrString += "\tCurrent Delta T         = " + fp2str(dt) + "\n";
    }
    if (unstablePhase) {
        ErrString += "\t  - The phase became unstable (Cp < 0) T_unstable_last = "
                     + fp2str(Tunstable) + "\n";
    }
    if (doUV) {
        throw CanteraError("setState_HPorUV (UV)", ErrString);
    } else {
        throw CanteraError("setState_HPorUV (HP)", ErrString);
    }
}

void ThermoPhase::setState_SP(doublereal Starget, doublereal p,
                              doublereal dTtol)
{
    setState_SPorSV(Starget, p, dTtol, false);
}

void ThermoPhase::setState_SV(doublereal Starget, doublereal v,
                              doublereal dTtol)
{
    setState_SPorSV(Starget, v, dTtol, true);
}

void ThermoPhase::setState_SPorSV(doublereal Starget, doublereal p,
                                  doublereal dTtol, bool doSV)
{
    doublereal v = 0.0;
    doublereal dt;
    if (doSV) {
        v = p;
        if (v < 1.0E-300) {
            throw CanteraError("setState_SPorSV (SV)",
                               "Input specific volume is too small or negative. v = " + fp2str(v));
        }
        setDensity(1.0/v);
    } else {
        if (p < 1.0E-300) {
            throw CanteraError("setState_SPorSV (SP)",
                               "Input pressure is too small or negative. p = " + fp2str(p));
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
    double Told = Tnew;
    double Sold = Snew;

    bool ignoreBounds = false;
    // Unstable phases are those for which
    // Cp < 0.0. These are possible for cases where
    // we have passed the spinodal curve.
    bool unstablePhase = false;
    double Tunstable = -1.0;
    bool unstablePhaseNew = false;

    // Newton iteration
    for (int n = 0; n < 500; n++) {
        Told = Tnew;
        Sold = Snew;
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
        double denom = std::max(fabs(Starget), acpd * dTtol);
        double SConvErr = fabs((Serr * Tnew)/denom);
        if (SConvErr < 0.00001 *dTtol || fabs(dt) < dTtol) {
            return;
        }
    }
    // We are here when there hasn't been convergence
    /*
     * Formulate a detailed error message, since questions seem to
     * arise often about the lack of convergence.
     */
    string ErrString =  "No convergence in 500 iterations\n";
    if (doSV) {
        ErrString += "\tTarget Entropy          = " + fp2str(Starget) + "\n";
        ErrString += "\tCurrent Specific Volume = " + fp2str(v) + "\n";
        ErrString += "\tStarting Temperature    = " + fp2str(Tinit) + "\n";
        ErrString += "\tCurrent Temperature     = " + fp2str(Tnew) + "\n";
        ErrString += "\tCurrent Entropy          = " + fp2str(Snew) + "\n";
        ErrString += "\tCurrent Delta T         = " + fp2str(dt) + "\n";
    } else {
        ErrString += "\tTarget Entropy          = " + fp2str(Starget) + "\n";
        ErrString += "\tCurrent Pressure        = " + fp2str(p) + "\n";
        ErrString += "\tStarting Temperature    = " + fp2str(Tinit) + "\n";
        ErrString += "\tCurrent Temperature     = " + fp2str(Tnew) + "\n";
        ErrString += "\tCurrent Entropy         = " + fp2str(Snew) + "\n";
        ErrString += "\tCurrent Delta T         = " + fp2str(dt) + "\n";
    }
    if (unstablePhase) {
        ErrString += "\t  - The phase became unstable (Cp < 0) T_unstable_last = "
                     + fp2str(Tunstable) + "\n";
    }
    if (doSV) {
        throw CanteraError("setState_SPorSV (SV)", ErrString);
    } else {
        throw CanteraError("setState_SPorSV (SP)", ErrString);
    }
}

doublereal ThermoPhase::err(const std::string& msg) const
{
    throw CanteraError("ThermoPhase","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0.0;
}

void ThermoPhase::getUnitsStandardConc(double* uA, int k, int sizeUA) const
{
    for (int i = 0; i < sizeUA; i++) {
        if (i == 0) {
            uA[0] = 1.0;
        }
        if (i == 1) {
            uA[1] = -int(nDim());
        }
        if (i == 2) {
            uA[2] = 0.0;
        }
        if (i == 3) {
            uA[3] = 0.0;
        }
        if (i == 4) {
            uA[4] = 0.0;
        }
        if (i == 5) {
            uA[5] = 0.0;
        }
    }
}

void ThermoPhase::setSpeciesThermo(SpeciesThermo* spthermo)
{
    if (m_spthermo) {
        if (m_spthermo != spthermo) {
            delete m_spthermo;
        }
    }
    m_spthermo = spthermo;
}

SpeciesThermo& ThermoPhase::speciesThermo(int k)
{
    if (!m_spthermo) {
        throw CanteraError("ThermoPhase::speciesThermo()",
                           "species reference state thermo manager was not set");
    }
    return *m_spthermo;
}

void ThermoPhase::initThermoFile(const std::string& inputFile,
                                 const std::string& id)
{

    if (inputFile.size() == 0) {
        throw CanteraError("ThermoPhase::initThermoFile",
                           "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("initThermoFile","could not open "
                           +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */
    //XML_Node& phaseNode_XML = xml();
    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    XML_Node* fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
        throw CanteraError("ThermoPhase::initThermo",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }
    //fxml_phase->copy(&phaseNode_XML);
    //initThermoXML(*fxml_phase, id);
    bool m_ok = importPhase(*fxml_phase, this);
    if (!m_ok) {
        throw CanteraError("ThermoPhase::initThermoFile","importPhase failed ");
    }
    delete fxml;
}

void ThermoPhase::initThermoXML(XML_Node& phaseNode, const std::string& id)
{

    /*
     * and sets the state
     */
    if (phaseNode.hasChild("state")) {
        XML_Node& stateNode = phaseNode.child("state");
        setStateFromXML(stateNode);
    }
    setReferenceComposition(0);
}

void ThermoPhase::setReferenceComposition(const doublereal* const x)
{
    xMol_Ref.resize(m_kk);
    if (x) {
        for (size_t k = 0; k < m_kk; k++) {
            xMol_Ref[k] = x[k];
        }
    } else {
        getMoleFractions(DATA_PTR(xMol_Ref));
    }
    double sum = -1.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += xMol_Ref[k];
    }
    if (fabs(sum) > 1.0E-11) {
        throw CanteraError("ThermoPhase::setReferenceComposition",
                           "input mole fractions don't sum to 1.0");
    }

}

void ThermoPhase::getReferenceComposition(doublereal* const x) const
{
    for (size_t k = 0; k < m_kk; k++) {
        x[k] = xMol_Ref[k];
    }
}

void ThermoPhase::initThermo()
{
    // Check to see that there is at least one species defined in the phase
    if (m_kk == 0) {
        throw CanteraError("ThermoPhase::initThermo()",
                           "Number of species is equal to zero");
    }
    xMol_Ref.resize(m_kk, 0.0);
}
void ThermoPhase::installSlavePhases(Cantera::XML_Node* phaseNode)
{
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

void ThermoPhase::setElementPotentials(const vector_fp& lambda)
{
    doublereal rrt = 1.0/(GasConstant* temperature());
    size_t mm = nElements();
    if (lambda.size() < mm) {
        throw CanteraError("setElementPotentials", "lambda too small");
    }
    if (!m_hasElementPotentials) {
        m_lambdaRRT.resize(mm);
    }
    for (size_t m = 0; m < mm; m++) {
        m_lambdaRRT[m] = lambda[m] * rrt;
    }
    m_hasElementPotentials = true;
}

bool ThermoPhase::getElementPotentials(doublereal* lambda) const
{
    doublereal rt = GasConstant* temperature();
    if (m_hasElementPotentials) {
        for (size_t m = 0; m < nElements(); m++) {
            lambda[m] =  m_lambdaRRT[m] * rt;
        }
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

    /*
     * Evaluate the current base activity coefficients if necessary
     */
    std::vector<double> ActCoeff_Base(m_kk);
    getActivityCoefficients(DATA_PTR(ActCoeff_Base));
    std::vector<double> Xmol_Base(m_kk);
    getMoleFractions(DATA_PTR(Xmol_Base));

    // Make copies of ActCoeff and Xmol_ for use in taking differences
    std::vector<double> ActCoeff(m_kk);
    std::vector<double> Xmol(m_kk);
    double v_totalMoles = 1.0;
    double TMoles_base = v_totalMoles;

    /*
     *  Loop over the columns species to be deltad
     */
    for (size_t j = 0; j < m_kk; j++) {
        /*
         * Calculate a value for the delta moles of species j
         * -> NOte Xmol_[] and Tmoles are always positive or zero
         *    quantities.
         * -> experience has shown that you always need to make the deltas greater than needed to
         *    change the other mole fractions in order to capture some effects.
         */
        double moles_j_base = v_totalMoles * Xmol_Base[j];
        deltaMoles_j = 1.0E-7 * moles_j_base + v_totalMoles * 1.0E-13 + 1.0E-150;
        /*
         * Now, update the total moles in the phase and all of the
         * mole fractions based on this.
         */
        v_totalMoles = TMoles_base + deltaMoles_j;
        for (size_t k = 0; k < m_kk; k++) {
            Xmol[k] = Xmol_Base[k] * TMoles_base / v_totalMoles;
        }
        Xmol[j] = (moles_j_base + deltaMoles_j) / v_totalMoles;

        /*
         * Go get new values for the activity coefficients.
         * -> Note this calls setState_PX();
         */
        setState_PX(pres, DATA_PTR(Xmol));
        getActivityCoefficients(DATA_PTR(ActCoeff));

        /*
         * Calculate the column of the matrix
         */
        double* const lnActCoeffCol = dlnActCoeffdlnN + ld * j;
        for (size_t k = 0; k < m_kk; k++) {
            lnActCoeffCol[k] = (2*moles_j_base + deltaMoles_j) *(ActCoeff[k] - ActCoeff_Base[k]) /
                               ((ActCoeff[k] + ActCoeff_Base[k]) * deltaMoles_j);
        }
        /*
         * Revert to the base case Xmol_, v_totalMoles
         */
        v_totalMoles = TMoles_base;
        Xmol = Xmol_Base;
    }
    /*
     * Go get base values for the activity coefficients.
     * -> Note this calls setState_TPX() again;
     * -> Just wanted to make sure that cantera is in sync
     *    with VolPhase after this call.
     */
    setState_PX(pres, DATA_PTR(Xmol_Base));
}

std::string ThermoPhase::report(bool show_thermo) const
{
    char p[800];
    string s = "";
    try {
        if (name() != "") {
            sprintf(p, " \n  %s:\n", name().c_str());
            s += p;
        }
        sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
        s += p;
        sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
        s += p;
        sprintf(p, "           density    %12.6g  kg/m^3\n", density());
        s += p;
        sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
        s += p;

        doublereal phi = electricPotential();
        if (phi != 0.0) {
            sprintf(p, "         potential    %12.6g  V\n", phi);
            s += p;
        }
        if (show_thermo) {
            sprintf(p, " \n");
            s += p;
            sprintf(p, "                          1 kg            1 kmol\n");
            s += p;
            sprintf(p, "                       -----------      ------------\n");
            s += p;
            sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n",
                    enthalpy_mass(), enthalpy_mole());
            s += p;
            sprintf(p, "   internal energy    %12.6g     %12.4g     J\n",
                    intEnergy_mass(), intEnergy_mole());
            s += p;
            sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n",
                    entropy_mass(), entropy_mole());
            s += p;
            sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n",
                    gibbs_mass(), gibbs_mole());
            s += p;
            sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n",
                    cp_mass(), cp_mole());
            s += p;
            try {
                sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n",
                        cv_mass(), cv_mole());
                s += p;
            } catch (CanteraError& err) {
                err.save();
                sprintf(p, " heat capacity c_v    <not implemented>       \n");
                s += p;
            }
        }

        size_t kk = nSpecies();
        vector_fp x(kk);
        vector_fp y(kk);
        vector_fp mu(kk);
        getMoleFractions(&x[0]);
        getMassFractions(&y[0]);
        getChemPotentials(&mu[0]);
        doublereal rt = GasConstant * temperature();
        //if (th.nSpecies() > 1) {

        if (show_thermo) {
            sprintf(p, " \n                           X     "
                    "            Y          Chem. Pot. / RT    \n");
            s += p;
            sprintf(p, "                     -------------     "
                    "------------     ------------\n");
            s += p;
            for (size_t k = 0; k < kk; k++) {
                if (x[k] > SmallNumber) {
                    sprintf(p, "%18s   %12.6g     %12.6g     %12.6g\n",
                            speciesName(k).c_str(), x[k], y[k], mu[k]/rt);
                } else {
                    sprintf(p, "%18s   %12.6g     %12.6g     \n",
                            speciesName(k).c_str(), x[k], y[k]);
                }
                s += p;
            }
        } else {
            sprintf(p, " \n                           X     "
                    "            Y\n");
            s += p;
            sprintf(p, "                     -------------"
                    "     ------------\n");
            s += p;
            for (size_t k = 0; k < kk; k++) {
                sprintf(p, "%18s   %12.6g     %12.6g\n",
                        speciesName(k).c_str(), x[k], y[k]);
                s += p;
            }
        }
    }
    //}
    catch (CanteraError& err) {
        err.save();
    }
    return s;
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
