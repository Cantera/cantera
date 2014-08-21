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
    operator=(right);
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




//==================================================================================================================
}
