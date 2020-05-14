/**
 * @file PlasmaKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAKINETICS_H
#define CT_PLASMAKINETICS_H

#include "GasKinetics.h"
#include "cantera/plasma/PlasmaPhase.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary plasma-phase chemistry.
 * @ingroup kinetics
 */
class PlasmaKinetics : public GasKinetics
{
public:
    //! @name Constructors and General Information
    //! @{

    //! Constructor.
    /*!
     *  @param thermo  Pointer to the gas ThermoPhase (optional)
     */
    PlasmaKinetics(thermo_t* thermo = 0);

    virtual std::string kineticsType() const {
        return "Plasma";
    }

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{
    virtual void init();
    //@}

    //! Update rates of the plasma reaction
    virtual void update_rates_T();

protected:
    virtual void addPlasmaReaction(PlasmaReaction& r);

    //! Reaction index of each plasma reaction
    std::vector<size_t> m_plasmaIndx;

    //! Process index of each plasma reaction;
    std::vector<size_t> m_plasmaProcessIndx;

    //! Pointer to the plasma phase
    PlasmaPhase* m_plasma;
};

}

#endif
