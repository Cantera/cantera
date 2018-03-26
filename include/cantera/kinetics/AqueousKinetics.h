/**
 * @file AqueousKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_AQUEOUSKINETICS_H
#define CT_AQUEOUSKINETICS_H

#include "BulkKinetics.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary aqueous-phase chemistry. This kinetics
 * manager implements standard mass-action reaction rate expressions for liquids
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * @deprecated To be removed after Cantera 2.4
 *
 * @ingroup kinetics
 */
class AqueousKinetics : public BulkKinetics
{
public:
    /// Constructor. Creates an empty reaction mechanism.
    AqueousKinetics(thermo_t* thermo = 0);

    virtual std::string kineticsType() const {
        return "Aqueous";
    }

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(doublereal* kfwd);
    void updateROP();

    //! Update temperature-dependent portions of reaction rates
    void _update_rates_T();

    //! Update properties that depend on concentrations.
    void _update_rates_C();

    //! Update the equilibrium constants in molar units.
    void updateKc();

    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
};
}

#endif
