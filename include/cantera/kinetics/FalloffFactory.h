/**
 *  @file FalloffFactory.h
 *  Parameterizations for reaction falloff functions. Used by classes
 *  that implement gas-phase kinetics (GasKinetics, GRI_30_Kinetics)
 *  (see \ref falloffGroup and class \link Cantera::Falloff Falloff\endlink).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_NEWFALLOFF_H
#define CT_NEWFALLOFF_H

#include "cantera/base/FactoryBase.h"
#include "cantera/base/ct_thread.h"
#include "cantera/base/smart_ptr.h"
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
class FalloffFactory : public FactoryBase
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static FalloffFactory* factory() {
        ScopedLock lock(falloff_mutex);
        if (!s_factory) {
            s_factory = new FalloffFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        ScopedLock lock(falloff_mutex);
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
     * @return    Returns a pointer to a new Falloff class.
     */
    virtual Falloff* newFalloff(int type, const vector_fp& c);

private:
    //! Pointer to the single instance of the factory
    static FalloffFactory* s_factory;

    //! default constructor, which is defined as private
    FalloffFactory() {}

    //!  Mutex for use when calling the factory
    static mutex_t falloff_mutex;
};

//! @copydoc FalloffFactory::newFalloff
shared_ptr<Falloff> newFalloff(int type, const vector_fp& c);

}
#endif
