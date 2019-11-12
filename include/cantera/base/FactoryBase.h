/**
 *  @file FactoryBase.h
 *  File contains the FactoryBase class declarations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FACTORY_BASE
#define CT_FACTORY_BASE

#include <vector>
#include <mutex>
#include <unordered_map>
#include <functional>
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

//! Base class for factories.
/*!
 * This class maintains a registry of all factories that derive from it, and
 * deletes them all when its static method deleteFactories is invoked.
 */
class FactoryBase
{
public:

    //! destructor
    virtual ~FactoryBase() {
    }

    //! static function that deletes all factories in the internal registry
    //! maintained in a static variable
    static void deleteFactories() {
        for (const auto& f : s_vFactoryRegistry) {
            f->deleteFactory();
        }
        s_vFactoryRegistry.clear();
    }

protected:

    //! Constructor.
    /*!
     * Adds the current object to the current static list
     */
    FactoryBase() {
        s_vFactoryRegistry.push_back(this);
    }

    //! Virtual abstract function that deletes the factory
    /*!
     *  This must be properly defined in child objects.
     */
    virtual void deleteFactory() = 0;

private:
    //! statically held list of Factories.
    static std::vector<FactoryBase*> s_vFactoryRegistry;
};

//! Factory class that supports registering functions to create objects
//!
//! Template arguments for the class are the base type created by the factory,
//! followed by the types of any arguments which need to be passed to the
//! functions used to create objects, e.g. arguments to the constructor.
template <class T, typename ... Args>
class Factory : public FactoryBase {
public:
    virtual ~Factory() {}

    //! Create an object using the object construction function corresponding to
    //! "name" and the provided constructor arguments
    T* create(const std::string& name, Args... args) {
        return m_creators.at(canonicalize(name))(args...);
    }

    //! Register a new object construction function
    void reg(const std::string& name, std::function<T*(Args...)> f) {
        m_creators[name] = f;
    }

    //! Add an alias for an existing registered type
    void addAlias(const std::string& original, const std::string& alias) {
        if (!m_creators.count(original)) {
            throw CanteraError("Factory::addAlias",
                "Name '{}' not registered", original);
        }
        m_synonyms[alias] = original;
    }

    //! Get the canonical name registered for a type
    std::string canonicalize(const std::string& name) {
        if (m_creators.count(name)) {
            return name;
        } else if (m_synonyms.count(name)) {
            return m_synonyms.at(name);
        } else if (m_deprecated_names.count(name)) {
            warn_deprecated(name,
                fmt::format("Use '{}' instead.", m_deprecated_names.at(name)));
            return m_deprecated_names.at(name);
        } else {
            throw CanteraError("Factory::canonicalize", "No such type: '{}'", name);
        }
    }

    //! Returns true if `name` is registered with this factory
    bool exists(const std::string& name) const {
        return m_creators.count(name) || m_synonyms.count(name);
    }

protected:
    //! Add a deprecated alias for an existing registered type
    void addDeprecatedAlias(const std::string& original,
                            const std::string& alias) {
        if (!m_creators.count(original)) {
            throw CanteraError("Factory::addDeprecatedAlias",
                "Name '{}' not registered", original);
        }
        m_deprecated_names[alias] = original;
    }

private:
    std::unordered_map<std::string, std::function<T*(Args...)>> m_creators;

    //! Map of synonyms to canonical names
    std::unordered_map<std::string, std::string> m_synonyms;

    //! Map of deprecated synonyms to canonical names. Use of these names will
    //! show a deprecation warning.
    std::unordered_map<std::string, std::string> m_deprecated_names;
};

}

#endif
