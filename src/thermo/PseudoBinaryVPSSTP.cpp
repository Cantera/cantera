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

#include <cstdio>

using namespace std;

namespace Cantera
{
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

ThermoPhase*
PseudoBinaryVPSSTP::duplMyselfAsThermoPhase() const
{
    return new PseudoBinaryVPSSTP(*this);
}

int PseudoBinaryVPSSTP::eosType() const
{
    return 0;
}

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

doublereal PseudoBinaryVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("PseudoBinaryVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}

void PseudoBinaryVPSSTP::initThermo()
{
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}

void  PseudoBinaryVPSSTP::initLengths()
{
    m_kk = nSpecies();
    moleFractions_.resize(m_kk);
}

void PseudoBinaryVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id_);
}

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
            } catch (CanteraError& e) {
                e.save();
                sprintf(p, " heat capacity c_v    <not implemented>       \n");
                s += p;
            }
        }

    } catch (CanteraError& e) {
        e.save();
    }
    return s;
}

}
