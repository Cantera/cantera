/**
 * @file ReactingVolDomain.h
 *
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef EXTRAGLOBALRXN_H
#define EXTRAGLOBALRXN_H

#include "cantera/kinetics/InterfaceKinetics.h"
#include <string>
#include <vector>

namespace Cantera
{

//!  Class describing an extra global reaction, which is defined as
//!  a linear combination of actuals reactions, global or mass-action, creating a global stoichiometric result
/*!
 *   This is useful for defining thermodynamics of global processes that occur
 *   on a surface or in a homogeneous phase.
 *
 *   The class is set up via the function setupElemRxnVector(RxnVector, specialSpecies) which defines
 *   the vector of stoichiometric coefficients representing the base reaction to combine in order to 
 *   achieve the global result that's to be calculated. specialSpecies is the index of the species
 *   within the kinetics object that is used to identify the global reaction. Rates of progress
 *   are defined in terms of the production rate of the special species.
 *   
 *   @deprecated Unfinished implementation to be removed after Cantera 2.2.
 */
class ExtraGlobalRxn
{

public:
    //!  Constructor takes a default kinetics pointer
    /*!
     *   @param[in]  k_ptr   Pointer to a Kinetics class that will be used as the basis
     *                       for constructing this class.
     */
    ExtraGlobalRxn(Kinetics* k_ptr);

    //! Destructor
    virtual ~ExtraGlobalRxn() {}

    void setupElemRxnVector(double* RxnVector,
                            int specialSpecies = -1);
    std::string reactionString();
    double deltaSpecValue(double* speciesVectorProperty);

    std::vector<int>& reactants();
    std::vector<int>& products();
    bool isReversible();

    double ROPValue(double* ROPKinVector);
    double FwdROPValue(double* FwdROPElemKinVector, double* RevROPElemKinVector);
    double RevROPValue(double* FwdROPElemKinVector, double* RevROPElemKinVector);

    double reactantStoichCoeff(int kKin);
    double productStoichCoeff(int kKin);
    bool m_ThisIsASurfaceRxn;
    double deltaRxnVecValue(double* rxnVectorProperty);

    //! This kinetics operator is associated with just one
    //! homogeneous phase, associated with tpList[0] phase
    /*!
     * Kinetics object pointer
     */
    Cantera::Kinetics* m_kinetics;

    //! This kinetics operator is associated with multiple
    //! homogeneous and surface phases.
    /*!
     * This object owns the Kinetics object
     */
    Cantera::InterfaceKinetics* m_InterfaceKinetics;

    int m_nKinSpecies;

    //! Number of reactants in the global reaction
    int m_nReactants;

    //!  Vector of reactants that make up the global reaction
    /*!
     *   This is a list of reactants using the kinetic species index
     */
    std::vector<int> m_Reactants;

    //!  Vector of reactant stoichiometries that make up the global reaction
    /*!
     *   This is a list of reactant stoichiometries. The species index is given in 
     *   the member m_Reactants using the kinetic species index.
     */
    std::vector<doublereal> m_ReactantStoich;

    int m_nProducts;
    std::vector<int> m_Products;
    std::vector<doublereal> m_ProductStoich;

    int m_nNetSpecies;
    std::vector<int> m_NetSpecies;
    std::vector<doublereal> m_netStoich;

    int m_nRxns;
    std::vector<doublereal> m_ElemRxnVector;

    int m_SpecialSpecies;
    bool m_SpecialSpeciesProduct;
    int m_SS_index;

    int iphaseKin;
    bool m_ok;
    bool m_reversible;


};
}
#endif
