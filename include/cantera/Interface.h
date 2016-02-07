/**
 * @file Interface.h
 *   Declaration and Definition for the class Interface.
 */
#ifndef CXX_INTERFACE
#define CXX_INTERFACE

#include "thermo.h"
#include "kinetics.h"

namespace Cantera
{
//! An interface between multiple bulk phases.
/*!
 * This class is defined mostly for convenience. It inherits both from SurfPhase
 * and InterfaceKinetics. It therefore represents a surface phase, and also acts
 * as the kinetics manager to manage reactions occurring on the surface,
 * possibly involving species from other phases.
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

    Interface(const Interface& ii) :
        SurfPhase(ii),
        InterfaceKinetics(ii),
        m_ok(ii.m_ok),
        m_r(ii.m_r) {
    }

    Interface& operator=(const Interface& right) {
        if (this == &right) {
            return *this;
        }
        SurfPhase::operator=(right);
        InterfaceKinetics::operator=(right);
        m_ok = right.m_ok;
        m_r = right.m_r;
        return *this;
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
 */
inline Interface* importInterface(const std::string& infile,
                                  const std::string& id,
                                  std::vector<ThermoPhase*> phases)
{
    return new Interface(infile, id, phases);
}

}

#endif
