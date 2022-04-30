//! @file Solution.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOLUTION_H
#define CT_SOLUTION_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class ThermoPhase;
class Kinetics;
class Transport;

//! A container class holding managers for all pieces defining a phase
class Solution
{
protected:
    Solution();

public:
    virtual ~Solution() {}
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
    virtual void setThermo(shared_ptr<ThermoPhase> thermo);

    //! Set the Kinetics object
    virtual void setKinetics(shared_ptr<Kinetics> kinetics);

    //! Set the Transport object
    virtual void setTransport(shared_ptr<Transport> transport);

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

    //! Add a phase adjacent to this phase. Usually this means a higher-dimensional
    //! phase that participates in reactions in this phase.
    void addAdjacent(shared_ptr<Solution> adjacent);

    //! Get the Solution object for an adjacent phase by index
    shared_ptr<Solution> adjacent(size_t i) {
        return m_adjacent.at(i);
    }

    //! Get the Solution object for an adjacent phase by name
    shared_ptr<Solution> adjacent(const std::string& name) {
        return m_adjacentByName.at(name);
    }

    //! Get the number of adjacent phases
    size_t nAdjacent() const {
         return m_adjacent.size();
    }

    AnyMap parameters(bool withInput=false) const;

    //! Access input data associated with header definition
    const AnyMap& header() const;
    AnyMap& header();

    //! Retrieve source used for object creation; usually an input file name
    const std::string source() const;

    //! Overwrite source (only required if object is not created using newSolution)
    void setSource(const std::string& source);

protected:
    shared_ptr<ThermoPhase> m_thermo;  //!< ThermoPhase manager
    shared_ptr<Kinetics> m_kinetics;  //!< Kinetics manager
    shared_ptr<Transport> m_transport;  //!< Transport manager

    //! Adjacent phases, for access by index
    std::vector<shared_ptr<Solution>> m_adjacent;

    //! Adjacent phases, for access by name
    std::map<std::string, shared_ptr<Solution>> m_adjacentByName;

    AnyMap m_header;  //!< Additional input fields; usually from a YAML input file
};

//! Create and initialize a new Solution from an input file
/*!
 * This constructor wraps newPhase(), newKinetics() and newTransportMgr() routines
 * for initialization.
 *
 * @param infile name of the input file
 * @param name   name of the phase in the file. If this is blank, the first phase
 *               in the file is used.
 * @param transport name of the transport model. If blank, the transport model specified
 *                  in the phase definition is used.
 * @param adjacent vector containing names of adjacent phases that participate in this
 *                 phases kinetics. If empty, adjacent phases will be instantiated based
 *                 on the phase definition.
 * @returns an initialized Solution object.
 */
shared_ptr<Solution> newSolution(const std::string& infile, const std::string& name,
    const std::string& transport, const std::vector<std::string>& adjacent);

//! Create and initialize a new Solution manager from an input file
/*!
 * This constructor wraps newPhase(), newKinetics() and
 * newTransportMgr() routines for initialization.
 *
 * @param infile name of the input file
 * @param name   name of the phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @param transport name of the transport model.
 * @param adjacent vector containing adjacent Solution objects. If empty, adjacent
 *                 phases will be instantiated based on the phase definition.
 * @returns an initialized Solution object.
 */
shared_ptr<Solution> newSolution(const std::string& infile,
                                 const std::string& name="",
                                 const std::string& transport="",
                                 const std::vector<shared_ptr<Solution>>& adjacent={});

//! Create and initialize a new Solution manager from AnyMap objects
/*!
 * This constructor wraps newPhase(), newKinetics() and
 * newTransportMgr() routines for initialization.
 *
 * @param phaseNode the node containing the phase definition (that is, thermo model,
 *     list of species, and initial state)
 * @param rootNode the root node of the tree containing the phase definition, which
 *     will be used as the default location from which to read species definitions.
 * @param transport name of the transport model.
 * @param adjacent vector containing adjacent Solution objects. If empty, adjacent
 *                 phases will be instantiated based on the phase definition.
 * @param related  vector of phases related to the same root Solution object. Used
 *                 internally by newSolution() when creating complex interfaces where
 *                 a phase may be adjacent to multiple other phases but should be
 *                 instantiated only once.
 * @returns an initialized Solution object.
 */
shared_ptr<Solution> newSolution(
    const AnyMap& phaseNode, const AnyMap& rootNode=AnyMap(),
    const std::string& transport="",
    const std::vector<shared_ptr<Solution>>& adjacent={},
    const std::map<std::string, shared_ptr<Solution>>& related={});

}

#endif
