/**
 *  @file RxnMolChange.cpp
 *
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef RXNMOLCHANGE_H
#define RXNMOLCHANGE_H


#include <vector>

namespace Cantera
{
class ExtraGlobalRxn;
class Kinetics;


//! Class that includes some bookeeping entries for a reaction or a global reaction defined on a surface
/*!
 *   Note that all indexes refer to a specific interfacial or homogeneous kinetics object. It does not
 *   refer to the Phase list indexes.
 *   @deprecated Unfinished implementation to be removed after Cantera 2.2.
 */
class RxnMolChange
{
public:
    //! Main constructor for the class
    /*!
     *  @param kinPtr   Pointer to the kinetics base class
     *  @param irxn     Specific reaction index.
     */
    RxnMolChange(Cantera::Kinetics* kinPtr, int irxn);

    //! Destructor
    ~RxnMolChange() {}

    //!  Constructor for the object if the object refers to a global reaction
    /*!
     *  @param kinPtr   Pointer to the kinetics base class
     *  @param egr      Specific reaction index.
     */
    RxnMolChange(Cantera::Kinetics* kinPtr, Cantera::ExtraGlobalRxn* egr);

    //! Vector of mole changes for each phase in the Kinetics object due to the current reaction
    /*!
     *  This is the sum of the product stoichiometric coefficient minum the reactant stoichioemtric coefficient
     *  for all the species in a phase.
     *  The index is over all the phases listed in the Kinetics object.
     */
    std::vector<double> m_phaseMoleChange;

    std::vector<double> m_phaseReactantMoles;

    std::vector<double> m_phaseProductMoles;

    //! Vector of mass changes for each phase in the Kinetics object due to the current reaction
    /*!
     *  This is the sum of the product stoichiometric coefficient minum the reactant stoichioemtric coefficient
     *  index multiplied by the molecular weight for all species in a phase.
     *  The index is over all of the phases listed in the Kinetics object.
     */
    std::vector<double> m_phaseMassChange;

    //! Vector of mass changes for each phase in the Kinetics object due to the current reaction
    /*!
     *  This is the sum of the product stoichiometric coefficient minum the reactant stoichioemtric coefficient
     *  index multiplied by the charg for all species in a phase.
     *  The index is over all of the phases listed in the Kinetics object.
     */
    std::vector<double> m_phaseChargeChange;

    //! Vector of phase types in the reaction
    /*!
     *   Collection of eosTypes for all phases in the kinetics object
     *  The index is over all of the phases listed in the Kinetics object.
     */
    std::vector<int>    m_phaseTypes;

    //! Vector of phase dimensions for the reaction
    /*!
     *  Collection of nDims for all phases in the kinetics object
     *  The index is over all of the phases listed in the Kinetics object.
     */
    std::vector<int>    m_phaseDims;

    //!  Number of phases in the kientics object
    int m_nPhases;

    //!  Shallow pointer pointing to the kinetics object
    Cantera::Kinetics* m_kinBase;

    //! Reaction number within the kinetics object
    /*!
     *  If this is neg 1, then this reaction refers to a global reaction
     *  specified by the m_egr pointer.
     */
    int m_iRxn;

    //! Maximum change in charge of any phase due to this reaction
    double m_ChargeTransferInRxn;

    //! Electrochemical beta parameter for the reaction
    double m_beta;

    //! Pointer to the specification of the global reaction
    /*!
     *  This is 0, if the class refers to a single reaction in the kinetics object
     */
    Cantera::ExtraGlobalRxn* m_egr;
};
}

#endif
