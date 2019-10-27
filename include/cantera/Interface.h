/**
 * @file Interface.h
 *   Declaration and Definition for the class Interface.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CXX_INTERFACE
#define CXX_INTERFACE

#pragma message("warning: Interface.h is deprecated and will be removed after Cantera 2.5.")

#include "thermo.h"
#include "kinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"

namespace Cantera
{
//! An interface between multiple bulk phases.
/*!
 * This class is defined mostly for convenience. It inherits both from SurfPhase
 * and InterfaceKinetics. It therefore represents a surface phase, and also acts
 * as the kinetics manager to manage reactions occurring on the surface,
 * possibly involving species from other phases.
 *
 * @deprecated To be removed after Cantera 2.5.
 *             Replaceable with Solution and/or SurfPhase/InterfaceKinetics.
 */
class Interface :
    public SurfPhase,
    public InterfaceKinetics
{
public:
    //! Constructor.
    /*!
     * Construct an Interface instance from a specification in an input file.
     *
     * @param infile  Cantera input file in CTI or CTML format.
     * @param id      Identification string to distinguish between multiple
     *     definitions within one input file.
     * @param otherPhases  Neighboring phases that may participate in the
     *     reactions on this interface. Don't include the surface phase
     */
    Interface(const std::string& infile, std::string id,
              std::vector<ThermoPhase*> otherPhases) :
        m_ok(false),
        m_r(0) {
        warn_deprecated("class Interface",
            "To be removed after Cantera 2.5. "
            "Replaceable with Solution and/or SurfPhase/InterfaceKinetics.");

        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }

        XML_Node* x = get_XML_Node("#"+id, m_r);
        if (!x) {
            throw CanteraError("Interface","error in get_XML_Node");
        }
        importPhase(*x, this);
        otherPhases.push_back(this);
        importKinetics(*x, otherPhases, this);
        m_ok = true;
    }

    //! Not operator
    bool operator!() {
        return !m_ok;
    }

    //! return whether the object has been instantiated
    bool ready() const {
        return m_ok;
    }

protected:
    //! Flag indicating that the object has been instantiated
    bool m_ok;

    //! XML_Node pointer to the XML File object that contains the Surface and
    //! the Interfacial Reaction object description
    XML_Node* m_r;
};

//! Import an instance of class Interface from a specification in an input file.
/*!
 *  This is the preferred method to create an Interface instance.
 *
 * @deprecated To be removed after Cantera 2.5. Replaceable with Solution.
 */
inline Interface* importInterface(const std::string& infile,
                                  const std::string& id,
                                  std::vector<ThermoPhase*> phases)
{
    warn_deprecated("importInterface", "To be removed after Cantera 2.5. "
                    "Replaceable with Solution.");
    return new Interface(infile, id, phases);
}

}

#endif
