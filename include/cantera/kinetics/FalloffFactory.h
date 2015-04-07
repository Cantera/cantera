/**
 *  @file FalloffFactory.h
 *  Parameterizations for reaction falloff functions. Used by classes
 *  that implement gas-phase kinetics (GasKinetics, GRI_30_Kinetics)
 *  (see \ref falloffGroup and class \link Cantera::Falloff Falloff\endlink).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_NEWFALLOFF_H
#define CT_NEWFALLOFF_H

#include "cantera/base/ct_defs.h"
#include "reaction_defs.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/base/ct_thread.h"

namespace Cantera
{

/**
 *     @defgroup falloffGroup  Falloff Parameterizations
 *    This section describes the parameterizations used
 *    to describe the fall-off in reaction rate constants
 *    due to intermolecular energy transfer.
 *
 *    @ingroup chemkinetics
 */

/**
 * Base class for falloff function calculators. Each instance of a subclass of
 * Falloff computes one falloff function. This base class implements the
 * trivial falloff function F = 1.0.
 *
 * @ingroup falloffGroup
 */
class Falloff
{
public:
    //! Default constructor is empty
    Falloff() {}

    //! default destructor is empty
    virtual ~Falloff() {}

    /**
     * Initialize. Must be called before any other method is invoked.
     *
     * @param c Vector of coefficients of the parameterization. The number and
     *     meaning of these coefficients is subclass-dependent.
     */
    virtual void init(const vector_fp& c) {}

    /**
     * Update the temperature-dependent portions of the falloff
     * function, if any. This method evaluates temperature-dependent
     * intermediate results and stores them in the 'work' array.
     * If not overloaded, the default behavior is to do nothing.
     * @param T Temperature [K].
     * @param work storage space for intermediate results.
     */
    virtual void updateTemp(doublereal T, doublereal* work) const {}

    /**
     * The falloff function. This is defined so that the
     * rate coefficient is
     *
     * \f[  k = F(Pr)\frac{Pr}{1 + Pr}. \f]
     *
     * Here \f$ Pr \f$ is the reduced pressure, defined by
     *
     * \f[
     * Pr = \frac{k_0 [M]}{k_\infty}.
     * \f]
     *
     * @param pr reduced pressure (dimensionless).
     * @param work array of size workSize() containing cached
     *             temperature-dependent intermediate results from a prior call
     *             to updateTemp.
     *
     * @return Returns the value of the falloff function \f$ F \f$ defined above
     */
    virtual doublereal F(doublereal pr, const doublereal* work) const {
        return 1.0;
    }

    //! The size of the work array required.
    virtual size_t workSize() {
        return 0;
    }
};

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
        ScopedLock lock(falloff_mutex) ;
        if (!s_factory) {
            s_factory = new FalloffFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        ScopedLock lock(falloff_mutex);
        if (s_factory) {
            delete s_factory;
            s_factory = 0;
        }
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
    static mutex_t falloff_mutex ;
};

}
#endif
