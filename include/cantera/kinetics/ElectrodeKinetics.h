/**
 * @file ElectrodeKinetics.h
 *
 * @ingroup chemkinetics
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_ELECTRODEKINETICS_H
#define CT_ELECTRODEKINETICS_H

#include "InterfaceKinetics.h"


namespace Cantera
{


//!  A kinetics manager for heterogeneous reaction mechanisms. The
//!  reactions are assumed to occur at a 2D interface between two 3D phases.
/*!
 *   This class is a slight addition to the InterfaceKinetics class, adding
 *   several concepts. First we explicitly identify the electrode and solution
 *   phases. We will also assume that there is an electron phase. 
 *
 *  @ingroup chemkinetics
 *  @deprecated Unfinished implementation to be removed after Cantera 2.2.
 */
class ElectrodeKinetics : public InterfaceKinetics
{
public:
    //! Constructor
    /*!
     * @param thermo The optional parameter may be used to initialize
     *               the object with one ThermoPhase object.
     *               HKM Note -> Since the interface kinetics
     *               object will probably require multiple ThermoPhase
     *               objects, this is probably not a good idea
     *               to have this parameter.
     */
    ElectrodeKinetics(thermo_t* thermo = 0);

    /// Destructor.
    virtual ~ElectrodeKinetics();

    //! Copy Constructor
    ElectrodeKinetics(const ElectrodeKinetics& right);

    //! Assignment operator
    ElectrodeKinetics& operator=(const ElectrodeKinetics& right);

    //! Duplication function
    /*!
     *  @param tpVector Vector of ThermoPhase pointers. These are shallow pointers to the
     *                  ThermoPhase objects that will comprise the phases for the new object.
     * 
     *   @return        Returns the duplicated object as the base class Kinetics object.
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    virtual int type() const;

    //!  Identify the metal phase and the electrons species
    /*!
     *   We fill in the internal variables, metalPhaseIndex_ and kElectronIndex_ here
     */
    void identifyMetalPhase();

    //! Internal routine that updates the Rates of Progress of the reactions
    /*!
     *  This is actually the guts of the functionality of the object
     */
    virtual void updateROP();

    virtual void determineFwdOrdersBV(ReactionData& rdata, std::vector<doublereal>& fwdFullorders);

    //void addGlobalReaction(ReactionData& r);

    double calcForwardROP_BV(size_t irxn, size_t iBeta, double ioc, double nStoich, double nu, doublereal ioNet);

    double calcForwardROP_BV_NoAct(size_t irxn, size_t iBeta,  double ioc, double nStoich, double nu, doublereal ioNet);

    bool getExchangeCurrentDensityFormulation(size_t irxn, doublereal& nStoich, doublereal& OCV, doublereal& io,
                                                doublereal& overPotential, doublereal& beta, doublereal& resistance);

    //! Calculate the open circuit voltage of a given reaction
    /*!
     *    If the reaction has no electron transport, then return 0.0
     *
     *  @param irxn   Reaction id
     */
    double openCircuitVoltage(size_t irxn);

    double calcCurrentDensity(double nu, double nStoich, double io, double beta, double temp, doublereal resistivity = 0.0) const;

    double solveCurrentRes(doublereal nu, doublereal nStoich, doublereal ioc, doublereal beta, doublereal temp,
	                   doublereal resistivity = 0.0, int iprob = 0) const;

    //! Prepare the class for the addition of reactions
    /*!
     *    (virtual from Kinetics)
     *    We determine the metal phase and solution phase here
     */ 
    virtual void init();

    virtual void finalize();

    //! Vector of additional information about each reaction
    /*!
     *  This vector contains information about the phase mole change for each reaction,
     *  for example.
     */
    std::vector<RxnMolChange*> rmcVector;

protected:
    //! Index of the metal phase in the list of phases for this kinetics object. This is the electron phase.
    size_t metalPhaseIndex_;

    //! Index of the solution phase in the list of phases for this surface
    size_t solnPhaseIndex_;

    //! Index of the electrons species in the list of species for this surface kinetics, if none set it to -1
    size_t kElectronIndex_;

    

};
}

#endif
