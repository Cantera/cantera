/**
 * @file ElectrodeKinetics.h
 *
 * @ingroup chemkinetics
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
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
 *   several concepts. First we explicity identify the electrode and solution
 *   phases. We will also assume that there is an electron phase. 
 *
 *  @ingroup chemkinetics
 */
class ElectrodeKinetics : public InterfaceKinetics
{
public:
    //! Constructor
    /*!
     * @param thermo The optional parameter may be used to initialize
     *               the object with one ThermoPhase object.
     *               HKM Note -> Since the interface kinetics
     *               object will probably require multiple thermophase
     *               objects, this is probably not a good idea
     *               to have this parameter.
     */
    ElectrodeKinetics(thermo_t* thermo = 0);

    /// Destructor.
    virtual ~ElectrodeKinetics();

    //! Copy Constructor for the %Kinetics object.
    ElectrodeKinetics(const ElectrodeKinetics& right);

    //! Assignment operator
    ElectrodeKinetics& operator=(const ElectrodeKinetics& right);

    //! Duplication function
    /*!
     *  @param tpVector Vector of %ThermoPhase pointers. These are shallow pointers to the
     *                  %ThermoPhase objects that will comprise the phases for the new object.
     * 
     *   @return        Returns the duplicated object as the base class %Kinetics object.
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    virtual int type() const;

    //!  Identify the metal phase and the electrons species
    /*!
     *   We fill in the internal variables, metalPhaseRS_ and kElectronRS_ here
     */
    void identifyMetalPhase();





    //void addGlobalReaction(ReactionData& r);





protected:

    //! index of the metal phase in the list of phases for this surface
    size_t metalPhaseRS_;

    size_t electronPhaseRS_;

    //! Index of the solution phase in the list of phases for this surface
    size_t solnPhaseRS_;

    //! Index of the electrons species in the list of species for this surface kinetics, if none set it to -1
    size_t kElectronRS_;




};
}

#endif
