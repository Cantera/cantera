/**
 *  @file MolalityVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ molality based activity coefficient formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MolalityVPSSTP MolalityVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon activities
 * based on the molality scale.  These include most of the methods for
 * calculating liquid electrolyte thermodynamics.
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/base/stringUtils.h"

#include <iomanip>
#include <cstdio>
#include <fstream>

using namespace std;

namespace Cantera
{

/*
 * Default constructor.
 *
 * This doesn't do much more than initialize constants with
 * default values for water at 25C. Water molecular weight
 * comes from the default elements.xml file. It actually
 * differs slightly from the IAPWS95 value of 18.015268. However,
 * density conservation and therefore element conservation
 * is the more important principle to follow.
 */
MolalityVPSSTP::MolalityVPSSTP() :
    VPStandardStateTP(),
    m_indexSolvent(0),
    m_pHScalingType(PHSCALE_PITZER),
    m_indexCLM(npos),
    m_weightSolvent(18.01528),
    m_xmolSolventMIN(0.01),
    m_Mnaught(18.01528E-3)
{
    /*
     * Change the default to be that charge neutrality in the
     * phase is necessary condition for the proper specification
     * of thermodynamic functions within the phase
     */
    m_chargeNeutralityNecessary = true;
}

/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor
 */
MolalityVPSSTP::MolalityVPSSTP(const MolalityVPSSTP& b) :
    VPStandardStateTP(),
    m_indexSolvent(b.m_indexSolvent),
    m_pHScalingType(b.m_pHScalingType),
    m_indexCLM(b.m_indexCLM),
    m_xmolSolventMIN(b.m_xmolSolventMIN),
    m_Mnaught(b.m_Mnaught),
    m_molalities(b.m_molalities)
{
    *this = operator=(b);
}

/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
MolalityVPSSTP& MolalityVPSSTP::
operator=(const MolalityVPSSTP& b)
{
    if (&b != this) {
        VPStandardStateTP::operator=(b);
        m_indexSolvent     = b.m_indexSolvent;
        m_pHScalingType    = b.m_pHScalingType;
        m_indexCLM         = b.m_indexCLM;
        m_weightSolvent    = b.m_weightSolvent;
        m_xmolSolventMIN   = b.m_xmolSolventMIN;
        m_Mnaught          = b.m_Mnaught;
        m_molalities       = b.m_molalities;
    }
    return *this;
}

/**
 *
 * ~MolalityVPSSTP():   (virtual)
 *
 * Destructor: does nothing:
 *
 */
MolalityVPSSTP::~MolalityVPSSTP()
{
}

/*
 * This routine duplicates the current object and returns
 * a pointer to ThermoPhase.
 */
ThermoPhase*
MolalityVPSSTP::duplMyselfAsThermoPhase() const
{
    return new MolalityVPSSTP(*this);
}

/*
 *  -------------- Utilities -------------------------------
 */

// Equation of state type flag.
/*
 * The ThermoPhase base class returns
 * zero. Subclasses should define this to return a unique
 * non-zero value. Known constants defined for this purpose are
 * listed in mix_defs.h. The MolalityVPSSTP class also returns
 * zero, as it is a non-complete class.
 */
int MolalityVPSSTP::eosType() const
{
    return 0;
}

// Set the pH scale, which determines the scale for single-ion activity
// coefficients.
/*
 *  Single ion activity coefficients are not unique in terms of the
 *  representing actual measurable quantities.
 */
void MolalityVPSSTP::setpHScale(const int pHscaleType)
{
    m_pHScalingType = pHscaleType;
    if (pHscaleType !=  PHSCALE_PITZER && pHscaleType !=  PHSCALE_NBS) {
        throw CanteraError("MolalityVPSSTP::setpHScale",
                           "Unknown scale type: " + int2str(pHscaleType));
    }
}

// Reports the pH scale, which determines the scale for single-ion activity
// coefficients.
/*
 *  Single ion activity coefficients are not unique in terms of the
 *  representing actual measurable quantities.
 */
int MolalityVPSSTP::pHScale() const
{
    return m_pHScalingType;
}

/*
 * setSolvent():
 *  Utilities for Solvent ID and Molality
 *  Here we also calculate and store the molecular weight
 *  of the solvent and the m_Mnaught parameter.
 *  @param k index of the solvent.
 */
void MolalityVPSSTP::setSolvent(size_t k)
{
    if (k >= m_kk) {
        throw CanteraError("MolalityVPSSTP::setSolute ",
                           "bad value");
    }
    m_indexSolvent = k;
    AssertThrowMsg(m_indexSolvent==0, "MolalityVPSSTP::setSolvent",
                   "Molality-based methods limit solvent id to being 0");
    m_weightSolvent = molecularWeight(k);
    m_Mnaught = m_weightSolvent / 1000.;
}

/*
 * return the solvent id index number.
 */
size_t MolalityVPSSTP::solventIndex() const
{
    return m_indexSolvent;
}

/*
 * Sets the minimum mole fraction in the molality formulation. The
 * minimum mole fraction must be in the range 0 to 0.9.
 */
void  MolalityVPSSTP::
setMoleFSolventMin(doublereal xmolSolventMIN)
{
    if (xmolSolventMIN <= 0.0) {
        throw CanteraError("MolalityVPSSTP::setSolute ", "trouble");
    } else if (xmolSolventMIN > 0.9) {
        throw CanteraError("MolalityVPSSTP::setSolute ", "trouble");
    }
    m_xmolSolventMIN = xmolSolventMIN;
}

/**
 * Returns the minimum mole fraction in the molality formulation.
 */
doublereal MolalityVPSSTP::moleFSolventMin() const
{
    return m_xmolSolventMIN;
}

/*
 * calcMolalities():
 *   We calculate the vector of molalities of the species
 *   in the phase and store the result internally:
 * \f[
 *     m_i = (n_i) / (1000 * M_o * n_{o,p})
 * \f]
 *    where
 *    - \f$ M_o \f$ is the molecular weight of the solvent
 *    - \f$ n_o \f$ is the mole fraction of the solvent
 *    - \f$ n_i \f$ is the mole fraction of the solute.
 *    - \f$ n_{o,p} = max (n_{o, min}, n_o) \f$
 *    - \f$ n_{o,min} \f$ = minimum mole fraction of solvent allowed
 *              in the denominator.
 */
void MolalityVPSSTP::calcMolalities() const
{
    getMoleFractions(DATA_PTR(m_molalities));
    double xmolSolvent = m_molalities[m_indexSolvent];
    if (xmolSolvent < m_xmolSolventMIN) {
        xmolSolvent = m_xmolSolventMIN;
    }
    double denomInv = 1.0/ (m_Mnaught * xmolSolvent);
    for (size_t k = 0; k < m_kk; k++) {
        m_molalities[k] *= denomInv;
    }
}

/*
 * getMolalities():
 *   We calculate the vector of molalities of the species
 *   in the phase
 * \f[
 *     m_i = (n_i) / (1000 * M_o * n_{o,p})
 * \f]
 *    where
 *    - \f$ M_o \f$ is the molecular weight of the solvent
 *    - \f$ n_o \f$ is the mole fraction of the solvent
 *    - \f$ n_i \f$ is the mole fraction of the solute.
 *    - \f$ n_{o,p} = max (n_{o, min}, n_o) \f$
 *    - \f$ n_{o,min} \f$ = minimum mole fraction of solvent allowed
 *              in the denominator.
 */
void MolalityVPSSTP::getMolalities(doublereal* const molal) const
{
    calcMolalities();
    for (size_t k = 0; k < m_kk; k++) {
        molal[k] = m_molalities[k];
    }
}

/*
 * setMolalities():
 *   We are supplied with the molalities of all of the
 *   solute species. We then calculate the mole fractions of all
 *   species and update the ThermoPhase object.
 *
 *     m_i = (n_i) / (W_o/1000 * n_o_p)
 *
 *    where M_o is the molecular weight of the solvent
 *    n_o is the mole fraction of the solvent
 *    n_i is the mole fraction of the solute.
 *    n_o_p = max (n_o_min, n_o)
 *    n_o_min = minimum mole fraction of solvent allowed
 *              in the denominator.
 */
void MolalityVPSSTP::setMolalities(const doublereal* const molal)
{

    double Lsum = 1.0 / m_Mnaught;
    for (size_t k = 1; k < m_kk; k++) {
        m_molalities[k] = molal[k];
        Lsum += molal[k];
    }
    double tmp = 1.0 / Lsum;
    m_molalities[m_indexSolvent] = tmp / m_Mnaught;
    double sum = m_molalities[m_indexSolvent];
    for (size_t k = 1; k < m_kk; k++) {
        m_molalities[k] = tmp * molal[k];
        sum += m_molalities[k];
    }
    if (sum != 1.0) {
        tmp = 1.0 / sum;
        for (size_t k = 0; k < m_kk; k++) {
            m_molalities[k] *= tmp;
        }
    }
    setMoleFractions(DATA_PTR(m_molalities));
    /*
     * Essentially we don't trust the input: We calculate
     * the molalities from the mole fractions that we
     * just obtained.
     */
    calcMolalities();
}

/*
 * setMolalitiesByName()
 *
 *  This routine sets the molalities by name
 *  HKM -> Might need to be more complicated here, setting
 *         neutrals so that the existing mole fractions are
 *         preserved.
 */
void MolalityVPSSTP::setMolalitiesByName(compositionMap& mMap)
{
    size_t kk = nSpecies();
    doublereal x;
    /*
     * Get a vector of mole fractions
     */
    vector_fp mf(kk, 0.0);
    getMoleFractions(DATA_PTR(mf));
    double xmolS = mf[m_indexSolvent];
    double xmolSmin = std::max(xmolS, m_xmolSolventMIN);
    compositionMap::iterator p;
    for (size_t k = 0; k < kk; k++) {
        p = mMap.find(speciesName(k));
        if (p != mMap.end()) {
            x = mMap[speciesName(k)];
            if (x > 0.0) {
                mf[k] = x * m_Mnaught * xmolSmin;
            }
        }
    }
    /*
     * check charge neutrality
     */
    size_t largePos = npos;
    double cPos = 0.0;
    size_t largeNeg = npos;
    double cNeg = 0.0;
    double sum = 0.0;
    for (size_t k = 0; k < kk; k++) {
        double ch = charge(k);
        if (mf[k] > 0.0) {
            if (ch > 0.0) {
                if (ch * mf[k] > cPos) {
                    largePos = k;
                    cPos = ch * mf[k];
                }
            }
            if (ch < 0.0) {
                if (fabs(ch) * mf[k] > cNeg) {
                    largeNeg = k;
                    cNeg = fabs(ch) * mf[k];
                }
            }
        }
        sum += mf[k] * ch;
    }
    if (sum != 0.0) {
        if (sum > 0.0) {
            if (cPos > sum) {
                mf[largePos] -= sum / charge(largePos);
            } else {
                throw CanteraError("MolalityVPSSTP:setMolalitiesbyName",
                                   "unbalanced charges");
            }
        } else {
            if (cNeg > (-sum)) {
                mf[largeNeg] -= (-sum) / fabs(charge(largeNeg));
            } else {
                throw CanteraError("MolalityVPSSTP:setMolalitiesbyName",
                                   "unbalanced charges");
            }
        }

    }
    sum = 0.0;
    for (size_t k = 0; k < kk; k++) {
        sum += mf[k];
    }
    sum = 1.0/sum;
    for (size_t k = 0; k < kk; k++) {
        mf[k] *= sum;
    }
    setMoleFractions(DATA_PTR(mf));
    /*
     * After we formally set the mole fractions, we
     * calculate the molalities again and store it in
     * this object.
     */
    calcMolalities();
}

/*
 * setMolalitiesByNames()
 *
 *   Set the molalities of the solutes by name
 */
void MolalityVPSSTP::setMolalitiesByName(const std::string& x)
{
    compositionMap xx = parseCompString(x, speciesNames());
    setMolalitiesByName(xx);
}


/*
 * ------------ Molar Thermodynamic Properties ----------------------
 */


/*
 * - Activities, Standard States, Activity Concentrations -----------
 */

/*
 * This method returns the activity convention.
 * Currently, there are two activity conventions
 *  Molar-based activities
 *       Unit activity of species at either a hypothetical pure
 *       solution of the species or at a hypothetical
 *       pure ideal solution at infinite dilution
 *   cAC_CONVENTION_MOLAR 0
 *      - default
 *
 *  Molality based activities
 *       (unit activity of solutes at a hypothetical 1 molal
 *        solution referenced to infinite dilution at all
 *        pressures and temperatures).
 *       (solvent is still on molar basis).
 *   cAC_CONVENTION_MOLALITY 1
 *
 *  We set the convention to molality here.
 */
int MolalityVPSSTP::activityConvention() const
{
    return cAC_CONVENTION_MOLALITY;
}

void MolalityVPSSTP::getActivityConcentrations(doublereal* c) const
{
    err("getActivityConcentrations");
}

doublereal MolalityVPSSTP::standardConcentration(size_t k) const
{
    err("standardConcentration");
    return -1.0;
}

doublereal MolalityVPSSTP::logStandardConc(size_t k) const
{
    err("logStandardConc");
    return -1.0;
}

void MolalityVPSSTP::getActivities(doublereal* ac) const
{
    err("getActivities");
}

/*
 * Get the array of non-dimensional activity coefficients at
 * the current solution temperature, pressure, and
 * solution concentration.
 * These are mole fraction based activity coefficients. In this
 * object, their calculation is based on translating the values
 * of Molality based activity coefficients.
 *  See Denbigh p. 278 for a thorough discussion.
 *
 * Note, the solvent is treated differently. getMolalityActivityCoeff()
 * returns the molar based solvent activity coefficient already.
 * Therefore, we do not have to divide by x_s here.
 */
void MolalityVPSSTP::getActivityCoefficients(doublereal* ac) const
{
    getMolalityActivityCoefficients(ac);
    AssertThrow(m_indexSolvent==0, "MolalityVPSSTP::getActivityCoefficients");
    double xmolSolvent = moleFraction(m_indexSolvent);
    if (xmolSolvent < m_xmolSolventMIN) {
        xmolSolvent = m_xmolSolventMIN;
    }
    for (size_t k = 1; k < m_kk; k++) {
        ac[k] /= xmolSolvent;
    }
}

// Get the array of non-dimensional molality based
//  activity coefficients at the current solution temperature,
//  pressure, and  solution concentration.
/*
 *  See Denbigh p. 278 for a thorough discussion. This class must be overwritten in
 *  classes which derive from %MolalityVPSSTP. This function takes over from the
 *  molar-based activity coefficient calculation, getActivityCoefficients(), in
 *  derived classes.
 *
 *  Note these activity coefficients have the current pH scale applied to them.
 *
 * @param acMolality Output vector containing the molality based activity coefficients.
 *                   length: m_kk.
 */
void MolalityVPSSTP::getMolalityActivityCoefficients(doublereal* acMolality) const
{
    getUnscaledMolalityActivityCoefficients(acMolality);
    applyphScale(acMolality);
}

/*
 * osmotic coefficient:
 *
 *  Calculate the osmotic coefficient of the solvent. Note there
 *  are lots of definitions of the osmotic coefficient floating
 *  around. We use the one defined in the Pitzer's book:
 *  (Activity Coeff in Electrolyte Solutions, K. S. Pitzer
 *   CRC Press, Boca Raton, 1991, p. 85, Eqn. 28).
 *
 *        Definition:
 *         - sum(m_i) * Mnaught * oc = ln(activity_solvent)
 */
doublereal MolalityVPSSTP::osmoticCoefficient() const
{
    /*
     * First, we calculate the activities all over again
     */
    vector_fp act(m_kk);
    getActivities(DATA_PTR(act));
    /*
     * Then, we calculate the sum of the solvent molalities
     */
    double sum = 0;
    for (size_t k = 1; k < m_kk; k++) {
        sum += std::max(m_molalities[k], 0.0);
    }
    double oc = 1.0;
    double lac = log(act[m_indexSolvent]);
    if (sum > 1.0E-200) {
        oc = - lac / (m_Mnaught * sum);
    }
    return oc;
}


void MolalityVPSSTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */


doublereal MolalityVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("MolalityVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}

/*
 * Returns the units of the standard and general concentrations
 * Note they have the same units, as their divisor is
 * defined to be equal to the activity of the kth species
 * in the solution, which is unitless.
 *
 * This routine is used in print out applications where the
 * units are needed. Usually, MKS units are assumed throughout
 * the program and in the XML input files.
 *
 * On return uA contains the powers of the units (MKS assumed)
 * of the standard concentrations and generalized concentrations
 * for the kth species.
 *
 *  uA[0] = kmol units - default  = 1
 *  uA[1] = m    units - default  = -nDim(), the number of spatial
 *                                dimensions in the Phase class.
 *  uA[2] = kg   units - default  = 0;
 *  uA[3] = Pa(pressure) units - default = 0;
 *  uA[4] = Temperature units - default = 0;
 *  uA[5] = time units - default = 0
 */
void MolalityVPSSTP::getUnitsStandardConc(double* uA, int k, int sizeUA) const
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

void MolalityVPSSTP::setToEquilState(const doublereal* lambda_RT)
{
    updateStandardStateThermo();
    err("setToEquilState");
}

/*
 * Set the thermodynamic state.
 */
void MolalityVPSSTP::setStateFromXML(const XML_Node& state)
{
    VPStandardStateTP::setStateFromXML(state);
    string comp = ctml::getChildValue(state,"soluteMolalities");
    if (comp != "") {
        setMolalitiesByName(comp);
    }
    if (state.hasChild("pressure")) {
        double p = ctml::getFloat(state, "pressure", "pressure");
        setPressure(p);
    }
}

/*
 * Set the temperature (K), pressure (Pa), and molalities
 * (gmol kg-1) of the solutes
 */
void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p,
                                  const doublereal* const molalities)
{
    setMolalities(molalities);
    setState_TP(t, p);
}

/*
 * Set the temperature (K), pressure (Pa), and molalities.
 */
void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, compositionMap& m)
{
    setMolalitiesByName(m);
    setState_TP(t, p);
}

/*
 * Set the temperature (K), pressure (Pa), and molality.
 */
void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, const std::string& m)
{
    setMolalitiesByName(m);
    setState_TP(t, p);
}


/*
 * @internal Initialize. This method is provided to allow
 * subclasses to perform any initialization required after all
 * species have been added. For example, it might be used to
 * resize internal work arrays that must have an entry for
 * each species.  The base class implementation does nothing,
 * and subclasses that do not require initialization do not
 * need to overload this method.  When importing a CTML phase
 * description, this method is called just prior to returning
 * from function importPhase.
 *
 * @see importCTML.cpp
 */
void MolalityVPSSTP::initThermo()
{
    initLengths();
    VPStandardStateTP::initThermo();

    /*
     * The solvent defaults to species 0
     */
    setSolvent(0);
    /*
     * Find the Cl- species
     */
    m_indexCLM = findCLMIndex();
}

//  Get the array of unscaled non-dimensional molality based
//  activity coefficients at the current solution temperature,
//  pressure, and  solution concentration.
/*
 *  See Denbigh p. 278 for a thorough discussion. This class must be overwritten in
 *  classes which derive from %MolalityVPSSTP. This function takes over from the
 *  molar-based activity coefficient calculation, getActivityCoefficients(), in
 *  derived classes.
 *
 * @param acMolality Output vector containing the molality based activity coefficients.
 *                   length: m_kk.
 */
void MolalityVPSSTP::getUnscaledMolalityActivityCoefficients(doublereal* acMolality) const
{
    err("getUnscaledMolalityActivityCoefficients");
}

//  Apply the current phScale to a set of activity Coefficients or activities
/*
 *  See the Eq3/6 Manual for a thorough discussion.
 *
 * @param acMolality input/Output vector containing the molality based
 *                   activity coefficients. length: m_kk.
 */
void MolalityVPSSTP::applyphScale(doublereal* acMolality) const
{
    err("applyphScale");
}

//  Returns the index of the Cl- species.
/*
 *  The Cl- species is special in the sense that its single ion
 *  molality-based activity coefficient is used in the specification
 *  of the pH scale for single ions. Therefore, we need to know
 *  what species index Cl- is. If the species isn't in the species
 *  list then this routine returns -1, and we can't use the NBS
 *  pH scale.
 *
 *  Right now we use a restrictive interpretation. The species
 *  must be named "Cl-". It must consist of exactly one Cl and one E
 *  atom.
 */
size_t MolalityVPSSTP::findCLMIndex() const
{
    size_t indexCLM = npos;
    size_t eCl = npos;
    size_t eE = npos;
    size_t ne = nElements();
    string sn;
    for (size_t e = 0; e < ne; e++) {
        sn = elementName(e);
        if (sn == "Cl" || sn == "CL") {
            eCl = e;
            break;
        }
    }
    // We have failed if we can't find the Cl element index
    if (eCl == npos) {
        return npos;
    }
    for (size_t e = 0; e < ne; e++) {
        sn = elementName(e);
        if (sn == "E" || sn == "e") {
            eE = e;
            break;
        }
    }
    // We have failed if we can't find the E element index
    if (eE == npos) {
        return npos;
    }
    for (size_t k = 1; k < m_kk; k++) {
        doublereal nCl = nAtoms(k, eCl);
        if (nCl != 1.0) {
            continue;
        }
        doublereal nE = nAtoms(k, eE);
        if (nE != 1.0) {
            continue;
        }
        for (size_t e = 0; e < ne; e++) {
            if (e != eE && e != eCl) {
                doublereal nA = nAtoms(k, e);
                if (nA != 0.0) {
                    continue;
                }
            }
        }
        sn = speciesName(k);
        if (sn != "Cl-" && sn != "CL-") {
            continue;
        }

        indexCLM = k;
        break;
    }
    return indexCLM;
}

//   Initialize lengths of local variables after all species have
//   been identified.
void  MolalityVPSSTP::initLengths()
{
    m_kk = nSpecies();
    m_molalities.resize(m_kk);
}

/*
 * initThermoXML()                (virtual from ThermoPhase)
 *   Import and initialize a ThermoPhase object
 *
 * @param phaseNode This object must be the phase node of a
 *             complete XML tree
 *             description of the phase, including all of the
 *             species data. In other words while "phase" must
 *             point to an XML phase object, it must have
 *             sibling nodes "speciesData" that describe
 *             the species in the phase.
 * @param id   ID of the phase. If nonnull, a check is done
 *             to see if phaseNode is pointing to the phase
 *             with the correct id.
 */
void MolalityVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{

    initLengths();
    /*
     * The solvent defaults to species 0
     */
    setSolvent(0);

    VPStandardStateTP::initThermoXML(phaseNode, id);
}

/**
  * Format a summary of the mixture state for output.
  */
std::string MolalityVPSSTP::report(bool show_thermo) const
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
        sprintf(p, "         potential    %12.6g  V\n", phi);
        s += p;

        size_t kk = nSpecies();
        vector_fp x(kk);
        vector_fp molal(kk);
        vector_fp mu(kk);
        vector_fp muss(kk);
        vector_fp acMolal(kk);
        vector_fp actMolal(kk);
        getMoleFractions(&x[0]);
        getMolalities(&molal[0]);
        getChemPotentials(&mu[0]);
        getStandardChemPotentials(&muss[0]);
        getMolalityActivityCoefficients(&acMolal[0]);
        getActivities(&actMolal[0]);

        size_t iHp = speciesIndex("H+");
        if (iHp != npos) {
            double pH = -log(actMolal[iHp]) / log(10.0);
            sprintf(p, "                pH    %12.4g  \n", pH);
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

        sprintf(p, " \n");
        s += p;
        if (show_thermo) {
            sprintf(p, "                           X        "
                    "   Molalities         Chem.Pot.    ChemPotSS    ActCoeffMolal\n");
            s += p;
            sprintf(p, "                                    "
                    "                      (J/kmol)      (J/kmol)                 \n");
            s += p;
            sprintf(p, "                     -------------  "
                    "  ------------     ------------  ------------    ------------\n");
            s += p;
            for (size_t k = 0; k < kk; k++) {
                if (x[k] > SmallNumber) {
                    sprintf(p, "%18s  %12.6g     %12.6g     %12.6g   %12.6g   %12.6g\n",
                            speciesName(k).c_str(), x[k], molal[k], mu[k], muss[k], acMolal[k]);
                } else {
                    sprintf(p, "%18s  %12.6g     %12.6g          N/A      %12.6g   %12.6g \n",
                            speciesName(k).c_str(), x[k], molal[k], muss[k], acMolal[k]);
                }
                s += p;
            }
        } else {
            sprintf(p, "                           X"
                    "Molalities\n");
            s += p;
            sprintf(p, "                     -------------"
                    "     ------------\n");
            s += p;
            for (size_t k = 0; k < kk; k++) {
                sprintf(p, "%18s   %12.6g     %12.6g\n",
                        speciesName(k).c_str(), x[k], molal[k]);
                s += p;
            }
        }
    } catch (CanteraError& err) {
        err.save();
    }
    return s;
}

/*
 * Format a summary of the mixture state for output.
 */
void MolalityVPSSTP::getCsvReportData(std::vector<std::string>& names,
                                      std::vector<vector_fp>& data) const
{
    names.clear();
    data.assign(10, vector_fp(nSpecies()));

    names.push_back("X");
    getMoleFractions(&data[0][0]);

    names.push_back("Molal");
    getMolalities(&data[1][0]);

    names.push_back("Chem. Pot. (J/kmol)");
    getChemPotentials(&data[2][0]);

    names.push_back("Chem. Pot. SS (J/kmol)");
    getStandardChemPotentials(&data[3][0]);

    names.push_back("Molal Act. Coeff.");
    getMolalityActivityCoefficients(&data[4][0]);

    names.push_back("Molal Activity");
    getActivities(&data[5][0]);

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
