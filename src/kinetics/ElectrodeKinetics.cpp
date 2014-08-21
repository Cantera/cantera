/**
 *  @file ElectrodeKinetics.cpp
 */

#include "cantera/kinetics/ElectrodeKinetics.h"

using namespace std;

namespace Cantera
{
//============================================================================================================================
ElectrodeKinetics::ElectrodeKinetics(thermo_t* thermo) :
    InterfaceKinetics(thermo),
    metalPhaseRS_(npos),
    electronPhaseRS_(npos),
    solnPhaseRS_(npos),
    kElectronRS_(npos)
{
 
}
//============================================================================================================================
ElectrodeKinetics::~ElectrodeKinetics()
{
  
}
//============================================================================================================================
ElectrodeKinetics::ElectrodeKinetics(const ElectrodeKinetics& right) :
    InterfaceKinetics()

{
    /*
     * Call the assignment operator
     */
    ElectrodeKinetics::operator=(right);
}
//============================================================================================================================
ElectrodeKinetics& ElectrodeKinetics::operator=(const ElectrodeKinetics& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    InterfaceKinetics::operator=(right);

    metalPhaseRS_   = right.metalPhaseRS_;
    electronPhaseRS_ = right.electronPhaseRS_;
    solnPhaseRS_ = right.solnPhaseRS_;
    kElectronRS_ = right.kElectronRS_;
   
    return *this;
}
//============================================================================================================================
int ElectrodeKinetics::type() const
{
    return cInterfaceKinetics;
}
//============================================================================================================================
Kinetics* ElectrodeKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    ElectrodeKinetics* iK = new ElectrodeKinetics(*this);
    iK->assignShallowPointers(tpVector);
    return iK;
}
//============================================================================================================================
//====================================================================================================================
//  Identify the metal phase and the electrons species
void ElectrodeKinetics::identifyMetalPhase()
{
    metalPhaseRS_ = npos;
    kElectronRS_ = -1;
    size_t np = nPhases();
    //
    // Identify the metal phase as the phase with the electron species (element index of 1 for element E
    // Should probably also stipulate a charge of -1.
    //
    for (size_t iph = 0; iph < np; iph++) {
        ThermoPhase* tp = m_thermo[iph];
        size_t nSpecies = tp->nSpecies();
        size_t nElements = tp->nElements();
        size_t eElectron = tp->elementIndex("E");
        if (eElectron != npos) {
            for (size_t k = 0; k < nSpecies; k++) {
                if (tp->nAtoms(k,eElectron) == 1) {
                    int ifound = 1;
                    for (size_t e = 0; e < nElements; e++) {
                        if (tp->nAtoms(k,e) != 0.0) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        metalPhaseRS_ = iph;
                        kElectronRS_ = m_start[iph] + k;
                    }
                }
            }
        }
       //
       //  Identify the solution phase as a 3D phase, with nonzero phase charge change 
       //  in at least one reaction
       //
       if (iph != metalPhaseRS_) {
            for (size_t i = 0; i < m_ii; i++) {
                RxnMolChange* rmc = rmcVector[i];
                if (rmc->m_phaseChargeChange[iph] != 0) {
                    if (rmc->m_phaseDims[iph] == 3) {
                       solnPhaseRS_ = iph;
                       break;
                    }
                }
             }
        }

    }
}

//==================================================================================================================
}
