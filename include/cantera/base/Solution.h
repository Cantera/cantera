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
class ExternalHandle;

//! @defgroup solnGroup Objects Representing Phases
//! High-level interface to %Cantera's core objects.

//! A container class for chemically-reacting solutions.
/*!
 * The Solution class collects all objects needed to describe a chemically-reacting
 * solution. Instances can be created to represent any type of solution -- a mixture
 * of gases, a liquid solution, or a solid solution, for example.
 *
 * Solution objects only define a small number of methods of their own, and are provided
 * so that a single object can be used to access thermodynamic, kinetic, and transport
 * properties of a solution:
 *  - ThermoPhase manager; accessed via thermo()
 *  - Kinetics manager; accessed via kinetics()
 *  - Transport manager; accessed via transport()
 *
 * The most common way to instantiate Solution objects is by using a phase definition,
 * species and reactions defined in an input file:
 * @code
 *    shared_ptr<Solution> sol = newSolution("gri30.yaml", "gri30");
 * @endcode
 * @ingroup solnGroup
 */
class Solution : public std::enable_shared_from_this<Solution>
{
protected:
    Solution() = default;

public:
    virtual ~Solution() = default;
    Solution(const Solution&) = delete;
    Solution& operator=(const Solution&) = delete;

    //! Create an empty Solution object
    static shared_ptr<Solution> create() {
        return shared_ptr<Solution>( new Solution );
    }

    //! Return the name of this Solution object
    string name() const;

    //! Set the name of this Solution object
    void setName(const string& name);

    //! Set the ThermoPhase object
    virtual void setThermo(shared_ptr<ThermoPhase> thermo);

    //! Set the Kinetics object
    virtual void setKinetics(shared_ptr<Kinetics> kinetics);

    //! Set the Transport object directly
    virtual void setTransport(shared_ptr<Transport> transport);

    //! Set the Transport object by name
    //! @param model  name of transport model; if omitted, the default model is used
    //! @since New in %Cantera 3.0
    void setTransportModel(const string& model="default");

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
    shared_ptr<Solution> adjacent(const string& name) {
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
    const string source() const;

    //! Overwrite source (only required if object is not created using newSolution)
    void setSource(const string& source);

    //! Store a handle to a wrapper for this Solution object from an external
    //! language interface (for example, a Python Solution object)
    void holdExternalHandle(const string& name, shared_ptr<ExternalHandle> handle);

    //! Get the handle for a wrapper for this Solution object from an external
    //! language interface.
    //! Returns a null pointer if the requested handle does not exist.
    shared_ptr<ExternalHandle> getExternalHandle(const string& name) const;

    //! Register a function to be called if any of the Solution's thermo, kinetics,
    //! or transport objects is replaced.
    //! @param id  A unique ID corresponding to the object affected by the callback.
    //!   Typically, this is a pointer to an object that also holds a reference to the
    //!   Solution object.
    //! @param callback  The callback function to be called after any component of the
    //!   Solution is replaced.
    //! When the callback becomes invalid (for example, the corresponding object is
    //! being deleted, the removeChangedCallback() method must be invoked.
    //! @since New in %Cantera 3.0
    void registerChangedCallback(void* id, const function<void()>& callback);

    //! Remove the callback function associated with the specified object.
    //! @since New in %Cantera 3.0
    void removeChangedCallback(void* id);

protected:
    shared_ptr<ThermoPhase> m_thermo;  //!< ThermoPhase manager
    shared_ptr<Kinetics> m_kinetics;  //!< Kinetics manager
    shared_ptr<Transport> m_transport;  //!< Transport manager

    //! Adjacent phases, for access by index
    vector<shared_ptr<Solution>> m_adjacent;

    //! Adjacent phases, for access by name
    map<string, shared_ptr<Solution>> m_adjacentByName;

    AnyMap m_header;  //!< Additional input fields; usually from a YAML input file

    //! Wrappers for this Solution object in extension languages, for evaluation
    //! of user-defined reaction rates
    map<string, shared_ptr<ExternalHandle>> m_externalHandles;

    //! Callback functions that are invoked when the therm, kinetics, or transport
    //! members of the Solution are replaced.
    map<void*, function<void()>> m_changeCallbacks;
};

//! Create and initialize a new Solution from an input file
/*!
 * This constructor wraps newThermo(), newKinetics() and newTransport() routines
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
 * @ingroup solnGroup
 */
shared_ptr<Solution> newSolution(const string& infile, const string& name,
    const string& transport, const vector<string>& adjacent);

//! Create and initialize a new Solution manager from an input file
/*!
 * This constructor wraps newThermo(), newKinetics() and newTransport() routines
 * for initialization.
 *
 * @param infile name of the input file
 * @param name   name of the phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @param transport name of the transport model.
 * @param adjacent vector containing adjacent Solution objects. If empty, adjacent
 *                 phases will be instantiated based on the phase definition.
 * @returns an initialized Solution object.
 * @ingroup solnGroup
 */
shared_ptr<Solution> newSolution(const string& infile,
                                 const string& name="",
                                 const string& transport="default",
                                 const vector<shared_ptr<Solution>>& adjacent={});

//! Create and initialize a new Solution manager from AnyMap objects
/*!
 * This constructor wraps newThermo(), newKinetics() and newTransport() routines
 * for initialization.
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
 * @ingroup solnGroup
 */
shared_ptr<Solution> newSolution(
    const AnyMap& phaseNode, const AnyMap& rootNode=AnyMap(),
    const string& transport="default",
    const vector<shared_ptr<Solution>>& adjacent={},
    const map<string, shared_ptr<Solution>>& related={});

}

#endif
