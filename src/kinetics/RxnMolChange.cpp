/**
 *  @file example2.cpp
 *
 * $Id: RxnMolChange.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 *
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/kinetics/RxnMolChange.h"


#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"

#include "cantera/kinetics/ExtraGlobalRxn.h"

#include <iostream>
#include <new>

using namespace Cantera;
using namespace std;

namespace Cantera {

RxnMolChange::RxnMolChange(Cantera::Kinetics* kinPtr, int irxn) :
    m_nPhases(0),
    m_kinBase(kinPtr),
    m_iRxn(irxn),
    m_ChargeTransferInRxn(0.0),
    m_beta(0.0),
    m_egr(0)
{
    warn_deprecated("class RxnMolChange", "To be removed after Cantera 2.2.");
    int iph;
    AssertTrace(irxn >= 0);
    AssertTrace(irxn < static_cast<int>(kinPtr->nReactions()));

    m_nPhases = static_cast<int>(kinPtr->nPhases());

    m_phaseMoleChange.resize(m_nPhases, 0.0);
    m_phaseReactantMoles.resize(m_nPhases, 0.0);
    m_phaseProductMoles.resize(m_nPhases, 0.0);
    m_phaseMassChange.resize(m_nPhases, 0.0);
    m_phaseChargeChange.resize(m_nPhases, 0.0);
    m_phaseTypes.resize(m_nPhases, 0);
    m_phaseDims.resize(m_nPhases, 0);

    int m_kk = static_cast<int>(kinPtr->nTotalSpecies());

    for (int kKin = 0; kKin < m_kk; kKin++) {
        iph = static_cast<int>(m_kinBase->speciesPhaseIndex(kKin));
        Cantera::ThermoPhase& tpRef = m_kinBase->thermo(iph);
        int kLoc = kKin - static_cast<int>(m_kinBase->kineticsSpeciesIndex(0, iph));
        double rsc = m_kinBase->reactantStoichCoeff(kKin, irxn);
        double psc = m_kinBase->productStoichCoeff(kKin, irxn);
        double nsc = psc - rsc;
        m_phaseMoleChange[iph] += (nsc);
        m_phaseReactantMoles[iph] += rsc;
        m_phaseProductMoles[iph] += psc;
        double mw = tpRef.molecularWeight(kLoc);
        m_phaseMassChange[iph] += (nsc) * mw;
        double chg = tpRef.charge(kLoc);
        m_phaseChargeChange[iph] += nsc * chg;
    }

    for (iph = 0; iph < m_nPhases; iph++) {
        Cantera::ThermoPhase& tpRef = m_kinBase->thermo(iph);
        m_phaseDims[iph] = static_cast<int>(tpRef.nDim());
        m_phaseTypes[iph] = tpRef.eosType();
        if (m_phaseChargeChange[iph] != 0.0) {
            double tmp = fabs(m_phaseChargeChange[iph]);
            if (tmp >  m_ChargeTransferInRxn) {
                m_ChargeTransferInRxn = tmp;
            }
        }
    }

    if (m_ChargeTransferInRxn) {
        Cantera::InterfaceKinetics* iK = dynamic_cast<Cantera::InterfaceKinetics*>(kinPtr);
        if (iK) {
            m_beta = iK->electrochem_beta(irxn);
        } else {
            throw Cantera::CanteraError("RxnMolChange", "unknown condition on charge");
        }
    }

}

RxnMolChange::RxnMolChange(Cantera::Kinetics* kinPtr, Cantera::ExtraGlobalRxn* egr) :
    m_nPhases(0),
    m_kinBase(kinPtr),
    m_iRxn(-1),
    m_ChargeTransferInRxn(0.0),
    m_beta(0.0),
    m_egr(egr)
{
    int iph;
    AssertTrace(egr != 0);

    m_nPhases = static_cast<int>(kinPtr->nPhases());

    m_phaseMoleChange.resize(m_nPhases, 0.0);
    m_phaseReactantMoles.resize(m_nPhases, 0.0);
    m_phaseProductMoles.resize(m_nPhases, 0.0);
    m_phaseMassChange.resize(m_nPhases, 0.0);
    m_phaseChargeChange.resize(m_nPhases, 0.0);
    m_phaseTypes.resize(m_nPhases, 0);
    m_phaseDims.resize(m_nPhases, 0);

    int m_kk = static_cast<int>(kinPtr->nTotalSpecies());

    for (int kKin = 0; kKin < m_kk; kKin++) {
        iph = static_cast<int>(m_kinBase->speciesPhaseIndex(kKin));
        ThermoPhase& tpRef = m_kinBase->thermo(iph);
        int kLoc = kKin - static_cast<int>(m_kinBase->kineticsSpeciesIndex(0, iph));
        double rsc = egr->reactantStoichCoeff(kKin);
        double psc = egr->productStoichCoeff(kKin);
        double nsc = psc - rsc;
        m_phaseMoleChange[iph] += (nsc);
        m_phaseReactantMoles[iph] += rsc;
        m_phaseProductMoles[iph] += psc;
        double mw = tpRef.molecularWeight(kLoc);
        m_phaseMassChange[iph] += (nsc) * mw;
        double chg = tpRef.charge(kLoc);
        m_phaseChargeChange[iph] += nsc * chg;
    }

    for (iph = 0; iph < m_nPhases; iph++) {
        ThermoPhase& tpRef = m_kinBase->thermo(iph);
        m_phaseDims[iph] = static_cast<int>(tpRef.nDim());
        m_phaseTypes[iph] = tpRef.eosType();
        if (m_phaseChargeChange[iph] != 0.0) {
            double tmp = fabs(m_phaseChargeChange[iph]);
            if (tmp >  m_ChargeTransferInRxn) {
                m_ChargeTransferInRxn = tmp;
            }
        }
    }

    if (m_ChargeTransferInRxn) {
        InterfaceKinetics* iK = dynamic_cast<InterfaceKinetics*>(kinPtr);
        if (iK) {
            m_beta = 0.0;
        } else {
            throw CanteraError("RxnMolChange", "unknown condition on charge");
        }
    }

}

}

