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
class ReactionFactory : public Factory<Reaction>
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

    void setup_XML(std::string name,
                   Reaction* R, const XML_Node& rxn_node) {
        try {
            m_xml_setup.at(name)(R, rxn_node);
        } catch (std::out_of_range&) {
            throw CanteraError("ReactionFactory::setup_XML",
                               "No such type: '{}'", name);
        }
    }

    //! Register a new object initializer function (from XML node)
    void reg_XML(const std::string& name,
                 std::function<void(Reaction*, const XML_Node&)> f) {
        m_xml_setup[name] = f;
    }

    void setup_AnyMap(std::string name,
                      Reaction* R, const AnyMap& node, const Kinetics& kin) {
        try {
            m_anymap_setup.at(name)(R, node, kin);
        } catch (std::out_of_range&) {
            throw CanteraError("ReactionFactory::setup_AnyMap",
                               "No such type: '{}'", name);
        }
    }

    //! Register a new object initializer function (from XML node)
    void reg_AnyMap(const std::string& name,
                    std::function<void(Reaction*, const AnyMap&, const Kinetics&)> f) {
        m_anymap_setup[name] = f;
    }

private:
    //! Pointer to the single instance of the factory
    static ReactionFactory* s_factory;

    //! default constructor, which is defined as private
    ReactionFactory();

    //!  Mutex for use when calling the factory
    static std::mutex reaction_mutex;

    //! map of XML initializers
    std::unordered_map<std::string,
        std::function<void(Reaction*, const XML_Node&)>> m_xml_setup;

    //! map of AnyMap initializers
    std::unordered_map<std::string,
        std::function<void(Reaction*, const AnyMap&, const Kinetics&)>> m_anymap_setup;
};

//! Create a new empty Reaction object
/*!
 * @param type string identifying type of reaction.
 */
unique_ptr<Reaction> newReaction(const std::string& type);

//! Create a new Reaction object for the reaction defined in `rxn_node`
/*!
 * @param rxn_node XML node describing reaction.
 */
unique_ptr<Reaction> newReaction(const XML_Node& rxn_node);

//! Create a new Reaction object using the specified parameters
/*!
 * @param node AnyMap node describing reaction.
 * @param kin kinetics manager
 */
unique_ptr<Reaction> newReaction(const AnyMap& node,
                                 const Kinetics& kin);
}
#endif
