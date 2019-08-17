//! @file Base.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BASE_H
#define CT_BASE_H

#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class ThermoPhase;
class Kinetics;
class Transport;

//! A container class holding managers for all pieces defining a phase
class SolutionBase : public std::enable_shared_from_this<SolutionBase>
{
public:
    SolutionBase();
    SolutionBase(const std::string& infile, const std::string& phasename);
    ~SolutionBase() {}
    SolutionBase(const SolutionBase&) = delete;
    SolutionBase& operator=(const SolutionBase&) = delete;

    static shared_ptr<SolutionBase> create() {
        return shared_ptr<SolutionBase>( new SolutionBase );
    }

    //! String indicating the type of the SolutionBase object. Corresponds
    //! to the type of phase originally instantiated
    std::string type() const;

    //! Set the type of this SolutionBase object
    void setType(const std::string& type);

    /*! Name and phase
     * Class SolutionBase references two strings that identify a SolutionBase.
     * The phase is the value of the phase name in YAML (or ID attribute
     * of the XML node) that is used to initialize a phase when it is read.
     * The name field is also initialized to the value of the phase name
     * read from the XML/YAML node.
     *
     * However, the name field may be changed to another value during the course
     * of a calculation. For example, if a SolutionBase is located in two places,
     * but has the same constitutive input, the phase of the two SolutionBases
     * will be the same, but the names of the two SolutionBases need to differ.
     *
     * It is an error to have two phases in a single problem with the same name
     * attribute. Thus, it is expected that there is a 1-1 correspondence between
     * names and unique SolutionBase objects within a Cantera problem.
     */

    //! Return the phase string of this SolutionBase.
    std::string phase() const;

    //! Set the phase string of this SolutionBase.
    void setPhase(const std::string& id);

    //! Return the name of this SolutionBase object
    std::string name() const;

    //! Set the name of this SolutionBase object
    void setName(const std::string& name);

    //! Generate self-documenting YAML string
    virtual std::string toYAML() const;

    //! Set the ThermoPhase object
    void setThermoPhase(shared_ptr<ThermoPhase> thermo);

    //! Set the Kinetics object
    void setKinetics(shared_ptr<Kinetics> kinetics);

    //! Set the Transport object
    void setTransport(shared_ptr<Transport> transport);

protected:
    shared_ptr<ThermoPhase> m_thermo;  //! ThermoPhase manager
    shared_ptr<Kinetics> m_kinetics;  //! Kinetics manager
    shared_ptr<Transport> m_transport;  //! Transport manager

    std::string m_type;  //! type of SolutionBase object
};

}
#endif
