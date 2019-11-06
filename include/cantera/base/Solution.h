//! @file Solution.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOLUTION_H
#define CT_SOLUTION_H

#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class ThermoPhase;
class Kinetics;
class Transport;

//! A container class holding managers for all pieces defining a phase
class Solution : public std::enable_shared_from_this<Solution>
{
private:
    Solution();

public:
    ~Solution() {}
    Solution(const Solution&) = delete;
    Solution& operator=(const Solution&) = delete;

    //! Create an empty Solution object
    static shared_ptr<Solution> create() {
        return shared_ptr<Solution>( new Solution );
    }

    //! Return the name of this Solution object
    std::string name() const;

    //! Set the name of this Solution object
    void setName(const std::string& name);

    //! Set the ThermoPhase object
    void setThermo(shared_ptr<ThermoPhase> thermo);

    //! Set the Kinetics object
    void setKinetics(shared_ptr<Kinetics> kinetics);

    //! Set the Transport object
    void setTransport(shared_ptr<Transport> transport);

    //! Accessor for the ThermoPhase pointer
    shared_ptr<ThermoPhase> thermo() {
        return m_thermo;
    }

    //! Accessor for the Kinetics pointer
    shared_ptr<Kinetics> kinetics() {
        return m_kinetics;
    }

    //! Accessor for the Transport pointer
    shared_ptr<Transport> transport() {
        return m_transport;
    }

protected:
    shared_ptr<ThermoPhase> m_thermo;  //!< ThermoPhase manager
    shared_ptr<Kinetics> m_kinetics;  //!< Kinetics manager
    shared_ptr<Transport> m_transport;  //!< Transport manager
};

//! Create and initialize a new Solution manager from an input file
/*!
 * This constructor wraps newPhase(), newKinetics() and
 * newTransportMgr() routines for initialization.
 *
 * @param infile name of the input file
 * @param name   name of the phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @param transport name of the transport model.
 * @param adjacent vector containing adjacent solution objects.
 * @returns an initialized Solution object.
 */
shared_ptr<Solution> newSolution(const std::string& infile,
                                 const std::string& name="",
                                 const std::string& transport="",
                                 const std::vector<shared_ptr<Solution>>& adjacent={});

}
#endif
