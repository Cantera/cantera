/**
 *  @file PseudoBinaryVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations
 *  (see \ref thermoprops
 * and class \link Cantera::PseudoBinaryVPSSTP PseudoBinaryVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon expressions
 * for the excess gibbs free energy expressed as a function of
 * the mole fractions.
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/thermo/PseudoBinaryVPSSTP.h"
#include "cantera/base/stringUtils.h"

#include <cmath>
#include <cstdio>

using namespace std;

namespace Cantera
{

/*
 * Default constructor.
 *
 */
PseudoBinaryVPSSTP::PseudoBinaryVPSSTP() :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0),
    cationPhase_(0),
    anionPhase_(0)
{
}

/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor
 */
PseudoBinaryVPSSTP::PseudoBinaryVPSSTP(const PseudoBinaryVPSSTP& b) :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0),
    cationPhase_(0),
    anionPhase_(0)
{
    *this = operator=(b);
}

/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
PseudoBinaryVPSSTP& PseudoBinaryVPSSTP::
operator=(const PseudoBinaryVPSSTP& b)
{
    if (&b != this) {
        GibbsExcessVPSSTP::operator=(b);
    }

    PBType_                     = b.PBType_;
    numPBSpecies_               = b.numPBSpecies_;
    indexSpecialSpecies_        = b.indexSpecialSpecies_;
    PBMoleFractions_            = b.PBMoleFractions_;
    cationList_                 = b.cationList_;
    numCationSpecies_           = b.numCationSpecies_;
    anionList_                  = b.anionList_;
    numAnionSpecies_            = b.numAnionSpecies_;
    passThroughList_            = b.passThroughList_;
    numPassThroughSpecies_      = b.numPassThroughSpecies_;
    neutralPBindexStart         = b.neutralPBindexStart;
    cationPhase_                = b.cationPhase_;
    anionPhase_                 = b.anionPhase_;
    moleFractionsTmp_           = b.moleFractionsTmp_;

    return *this;
}

/**
 *
 * ~PseudoBinaryVPSSTP():   (virtual)
 *
 * Destructor: does nothing:
 *
 */
PseudoBinaryVPSSTP::~PseudoBinaryVPSSTP()
{
}

/*
 * This routine duplicates the current object and returns
 * a pointer to ThermoPhase.
 */
thermo_t*
PseudoBinaryVPSSTP::duplMyselfAsThermoPhase() const
{
    return new PseudoBinaryVPSSTP(*this);
}

/*
 *  -------------- Utilities -------------------------------
 */


// Equation of state type flag.
/*
 * The ThermoPhase base class returns
 * zero. Subclasses should define this to return a unique
 * non-zero value. Known constants defined for this purpose are
 * listed in mix_defs.h. The PseudoBinaryVPSSTP class also returns
 * zero, as it is a non-complete class.
 */
int PseudoBinaryVPSSTP::eosType() const
{
    return 0;
}



/*
 * ------------ Molar Thermodynamic Properties ----------------------
 */


/*
 * - Activities, Standard States, Activity Concentrations -----------
 */


doublereal PseudoBinaryVPSSTP::standardConcentration(size_t k) const
{
    err("standardConcentration");
    return -1.0;
}

doublereal PseudoBinaryVPSSTP::logStandardConc(size_t k) const
{
    err("logStandardConc");
    return -1.0;
}



void PseudoBinaryVPSSTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}

void PseudoBinaryVPSSTP::calcPseudoBinaryMoleFractions() const
{
    size_t k;
    doublereal sumCat;
    doublereal sumAnion;
    doublereal sum = 0.0;
    switch (PBType_) {
    case PBTYPE_PASSTHROUGH:
        for (k = 0; k < m_kk; k++) {
            PBMoleFractions_[k] = moleFractions_[k];
        }
        break;
    case PBTYPE_SINGLEANION:
        sumCat = 0.0;
        sumAnion = 0.0;
        for (k = 0; k < m_kk; k++) {
            moleFractionsTmp_[k] = moleFractions_[k];
        }
        for (k = 0; k < cationList_.size(); k++) {
            sumCat += moleFractions_[cationList_[k]];
        }
        sumAnion =  moleFractions_[anionList_[k]];
        PBMoleFractions_[0] = sumCat -sumAnion;
        moleFractionsTmp_[indexSpecialSpecies_] -= PBMoleFractions_[0];


        for (k = 0; k < numCationSpecies_; k++) {
            PBMoleFractions_[1+k] = moleFractionsTmp_[cationList_[k]];
        }

        for (k = 0; k <  numPassThroughSpecies_; k++) {
            PBMoleFractions_[neutralPBindexStart + k] =
                moleFractions_[cationList_[k]];
        }

        sum = std::max(0.0, PBMoleFractions_[0]);
        for (k = 1; k < numPBSpecies_; k++) {
            sum += PBMoleFractions_[k];
        }
        for (k = 0; k < numPBSpecies_; k++) {
            PBMoleFractions_[k] /= sum;
        }

        break;
    case PBTYPE_SINGLECATION:
        throw CanteraError("eosType", "Unknown type");

        break;

    case PBTYPE_MULTICATIONANION:
        throw CanteraError("eosType", "Unknown type");

        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;

    }
}

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */


doublereal PseudoBinaryVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("PseudoBinaryVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
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
void PseudoBinaryVPSSTP::initThermo()
{
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}


//   Initialize lengths of local variables after all species have
//   been identified.
void  PseudoBinaryVPSSTP::initLengths()
{
    m_kk = nSpecies();
    moleFractions_.resize(m_kk);
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
void PseudoBinaryVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{


    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);
}

/**
  * Format a summary of the mixture state for output.
  */
std::string PseudoBinaryVPSSTP::report(bool show_thermo) const
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

        getChemPotentials(&mu[0]);
        getStandardChemPotentials(&muss[0]);
        getActivities(&actMolal[0]);


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

    } catch (CanteraError& err) {
        err.save();
    }
    return s;
}


}

