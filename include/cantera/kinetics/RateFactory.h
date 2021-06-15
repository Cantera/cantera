/**
 *  @file RateFactory.h
 *  Factory class for reaction rate objects. Used by classes
 *  that implement kinetics
 *  (see \ref reactionGroup and class \link Cantera::Rate Rate\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NEWRATE_H
#define CT_NEWRATE_H

#include "cantera/base/FactoryBase.h"
#include "cantera/kinetics/ReactionRate.h"

namespace Cantera
{

class Kinetics;
class Units;

/**
 * Factory class to construct reaction rate calculators.
 * The reaction factory is accessed through the static method factory:
 *
 * @code
 * Rate* f = RateFactory::factory()->newRate(type, c)
 * @endcode
 *
 * @ingroup reactionGroup
 */
class RateFactory : public Factory<ReactionRateBase, const AnyMap&, const Units&>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static RateFactory* factory() {
        std::unique_lock<std::mutex> lock(rate_mutex);
        if (!s_factory) {
            s_factory = new RateFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(rate_mutex);
        delete s_factory;
        s_factory = 0;
    }

private:
    //! Pointer to the single instance of the factory
    static RateFactory* s_factory;

    //! default constructor, which is defined as private
    RateFactory();

    //!  Mutex for use when calling the factory
    static std::mutex rate_mutex;
};


//! Create a new empty ReactionRateBase object
/*!
 * @param type string identifying type of reaction rate.
 */
shared_ptr<ReactionRateBase> newRate(const std::string& type);

//! Create a new Rate object using the specified parameters
/*!
 * @param rate_node AnyMap node describing reaction rate.
 * @param rate_units Unit system of the reaction rate
 */
shared_ptr<ReactionRateBase> newRate(const AnyMap& rate_node, const Units& rate_units);

//! Create a new Rate object using the specified parameters
/*!
 * @param rate_node AnyMap node describing reaction rate.
 * @param kin kinetics manager
 */
shared_ptr<ReactionRateBase> newRate(const AnyMap& rate_node, const Kinetics& kin);

}
#endif
