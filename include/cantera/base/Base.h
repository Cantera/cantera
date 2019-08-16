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
class SolutionBase {
public:
    SolutionBase();
    SolutionBase(const std::string& infile, const std::string& phasename);
    ~SolutionBase() {}
    SolutionBase(const SolutionBase&) = delete;
    SolutionBase& operator=(const SolutionBase&) = delete;

    //! String indicating the type of the SolutionBase object. Corresponds
    //! to the type of phase originally instantiated
    std::string type() const {
        return m_type;
    }

    //! Set the type of this SolutionBase object
    void setType(const std::string& type){
        m_type = type;
    }

    //! Return the name of this SolutionBase object
    std::string name() const {
        return m_name;
    }

    //! Set the name of this SolutionBase object
    void setName(const std::string& name){
        m_name = name;
    }

    //! Generate self-documenting YAML string
    virtual std::string toYAML() const {
        throw NotImplementedError("SolutionBase::toYAML");
    }

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
    std::string m_name;  //! name of SolutionBase object
};

}
#endif
