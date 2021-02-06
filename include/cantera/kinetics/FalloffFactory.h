/**
 *  @file FalloffFactory.h
 *  Parameterizations for reaction falloff functions. Used by classes
 *  that implement gas-phase kinetics (GasKinetics, GRI_30_Kinetics)
 *  (see \ref falloffGroup and class \link Cantera::Falloff Falloff\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NEWFALLOFF_H
#define CT_NEWFALLOFF_H

#include "cantera/base/FactoryBase.h"
#include "cantera/kinetics/Falloff.h"

namespace Cantera
{

/**
 * Factory class to construct falloff function calculators.
 * The falloff factory is accessed through static method factory:
 *
 * @code
 * Falloff* f = FalloffFactory::factory()->newFalloff(type, c)
 * @endcode
 *
 * @ingroup falloffGroup
 */
class FalloffFactory : public Factory<Falloff>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static FalloffFactory* factory() {
        std::unique_lock<std::mutex> lock(falloff_mutex);
        if (!s_factory) {
            s_factory = new FalloffFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(falloff_mutex);
        delete s_factory;
        s_factory = 0;
    }

    //! Return a pointer to a new falloff function calculator.
    /*!
     * @param type Integer flag specifying the type of falloff function. The
     *              standard types are defined in file reaction_defs.h. A
     *              factory class derived from FalloffFactory may define other
     *              types as well.
     * @param c    input vector of doubles which populates the falloff
     *             parameterization.
     * @returns    a pointer to a new Falloff class.
     *
     * @deprecated To be removed after Cantera 2.5.
     */
    virtual Falloff* newFalloff(int type, const vector_fp& c);

    //! Return a pointer to a new falloff function calculator.
    /*!
     * @param type String identifier specifying the type of falloff function.
     *             The standard types match class names defined in Falloff.h.
     *             A factory class derived from FalloffFactory may define
     *             other types as well.
     * @param c    input vector of doubles which populates the falloff
     *             parameterization.
     * @returns    a pointer to a new Falloff class.
     */
    virtual Falloff* newFalloff(const std::string& type, const vector_fp& c);

private:
    //! Pointer to the single instance of the factory
    static FalloffFactory* s_factory;

    //! default constructor, which is defined as private
    FalloffFactory();

    //!  Mutex for use when calling the factory
    static std::mutex falloff_mutex;
};

//! @copydoc FalloffFactory::newFalloff
/*!
 * @deprecated To be removed after Cantera 2.5.
 */
shared_ptr<Falloff> newFalloff(int type, const vector_fp& c);

//! @copydoc FalloffFactory::newFalloff
shared_ptr<Falloff> newFalloff(const std::string& type, const vector_fp& c);

}
#endif
