/**
 *  @file ReactionFactory.h
 *  Factory class for reaction functions. Used by classes
 *  that implement kinetics
 *  (see \ref reactionGroup and class \link Cantera::Reaction Reaction\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_NEWREACTION_H
#define CT_NEWREACTION_H

#include "cantera/base/FactoryBase.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

/**
 * Factory class to construct reaction function calculators.
 * The reaction factory is accessed through static method factory:
 *
 * @code
 * Reaction* f = ReactionFactory::factory()->newReaction(type, c)
 * @endcode
 *
 * @ingroup reactionGroup
 */
class ReactionFactory : public Factory<Reaction, const AnyMap&, const Kinetics&>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static ReactionFactory* factory() {
        std::unique_lock<std::mutex> lock(reaction_mutex);
        if (!s_factory) {
            s_factory = new ReactionFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(reaction_mutex);
        delete s_factory;
        s_factory = 0;
    }

private:
    //! Pointer to the single instance of the factory
    static ReactionFactory* s_factory;

    //! default constructor, which is defined as private
    ReactionFactory();

    //!  Mutex for use when calling the factory
    static std::mutex reaction_mutex;
};

class ReactionFactoryXML : public Factory<Reaction, const XML_Node&>
{
public:
    /**
     * Return a pointer to the factory. On the first call, a new instance is
     * created. Since there is no need to instantiate more than one factory,
     * on all subsequent calls, a pointer to the existing factory is returned.
     */
    static ReactionFactoryXML* factory() {
        std::unique_lock<std::mutex> lock(reaction_mutex);
        if (!s_factory) {
            s_factory = new ReactionFactoryXML;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(reaction_mutex);
        delete s_factory;
        s_factory = 0;
    }

private:
    //! Pointer to the single instance of the factory
    static ReactionFactoryXML* s_factory;

    //! default constructor, which is defined as private
    ReactionFactoryXML();

    //!  Mutex for use when calling the factory
    static std::mutex reaction_mutex;
};

}
#endif
