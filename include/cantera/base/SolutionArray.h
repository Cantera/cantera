//! @file SolutionArray.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOLUTIONARRAY_H
#define CT_SOLUTIONARRAY_H

#include "cantera/base/global.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class Solution;
class ThermoPhase;

/*!
 *  A container class providing a convenient interface for representing many
 *  thermodynamic states using the same Solution object. C++ SolutionArray objects are
 *  one-dimensional by design; extensions to multi-dimensional arrays need to be
 *  implemented in high-level API's.
 *
 *  @since  New in Cantera 3.0.
 *  @warning This function is an experimental part of the %Cantera API and may be
 *      changed or removed without notice.
 */
class SolutionArray
{
private:
    SolutionArray(const shared_ptr<Solution>& sol,
                  size_t size,
                  const AnyMap& meta);

public:
    virtual ~SolutionArray() {}

    /*!
     *  Instantiate a new SolutionArray reference
     *
     *  @param sol  Solution object defining phase definitions
     *  @param size  Number of SolutionArray entries
     *  @param meta  AnyMap holding SolutionArray meta data
     */
    static shared_ptr<SolutionArray> create(const shared_ptr<Solution>& sol,
                                            size_t size=0,
                                            const AnyMap& meta={})
    {
        return shared_ptr<SolutionArray>(new SolutionArray(sol, size, meta));
    }

    //! Size of SolutionArray (number of entries)
    int size() const {
        return m_size;
    }

    //! SolutionArray meta data.
    AnyMap& meta() {
        return m_meta;
    }

    //! Retrieve associated ThermoPhase object
    shared_ptr<ThermoPhase> thermo();

    /*!
     *  Check whether SolutionArray contains a component (property defining state or
     *  auxiliary variable)
     */
    bool hasComponent(const std::string& name) const;

    //! Retrieve a component of the SolutionArray by name
    vector_fp getComponent(const std::string& name) const;

    /*!
     *  Set a component of the SolutionArray by name.
     *
     *  @param name  Name of component (property defining state or auxiliary variable)
     *  @param data  Component data
     *  @param force  If true, add new component to SolutionArray
     */
    void setComponent(const std::string& name, const vector_fp& data, bool force=false);

    /*!
     *  Update the buffered index used to access entries.
     */
    void setIndex(size_t index, bool restore=true);

    /*!
     *  Retrieve the state vector for a single entry. If index is valid, it is updated;
     *  otherwise, the last previously used index is referenced.
     */
    vector_fp getState(size_t index=npos);

    /*!
     *  Set the state vector for a single entry. If index is valid, it is updated;
     *  otherwise, the last previously used index is referenced.
     */
    void setState(const vector_fp& data, size_t index=npos);

    /*!
     *  Retrieve auxiliary data for a single entry. If index is valid, it is updated;
     *  otherwise, the last previously used index is referenced.
     */
    std::map<std::string, double> getAuxiliary(size_t index=npos);

    /*!
     *  Write header data to container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of SolutionArray within the container file
     *  @param desc  Description
     */
    static void writeHeader(const std::string& fname, const std::string& id,
                            const std::string& desc);

    /*!
     *  Write header data to AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of SolutionArray node within AnyMap structure
     *  @param desc  Description
     */
    static void writeHeader(AnyMap& root, const std::string& id,
                            const std::string& desc);

    /*!
     *  Write SolutionArray data to container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of SolutionArray within the container file
     *  @param compression  Compression level; optional (default=0; HDF only)
     */
    void writeEntry(const std::string& fname, const std::string& id,
                    int compression=0);

    /*!
     *  Write SolutionArray data to AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of SolutionArray node within AnyMap structure
     */
    void writeEntry(AnyMap& root, const std::string& id);

    /*!
     *  Save current SolutionArray and header to a container file.
     *
     *  @param fname  Name of output container file (YAML or HDF)
     *  @param id  Identifier of SolutionArray within the container file
     *  @param desc  Description
     *  @param compression  Compression level; optional (default=0; HDF only)
     */
    void save(const std::string& fname, const std::string& id,
              const std::string& desc, int compression=0);

    /*!
     *  Read header data from container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of SolutionArray within the file structure
     */
    static AnyMap readHeader(const std::string& fname, const std::string& id);

    /*!
     *  Read header data from AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of SolutionArray node within AnyMap structure
     */
    static AnyMap readHeader(const AnyMap& root, const std::string& id);

    /*!
     *  Restore SolutionArray entry from a container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of SolutionArray within the file structure
     */
    void readEntry(const std::string& fname, const std::string& id);

    /*!
     *  Restore SolutionArray entry from AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of SolutionArray node within AnyMap structure
     */
    void readEntry(const AnyMap& root, const std::string& id);

    /*!
     *  Restore SolutionArray entry and header from a container file.
     *
     *  @param fname  Name of container file (YAML or HDF)
     *  @param id  Identifier of SolutionArray within the container file
     */
    AnyMap restore(const std::string& fname, const std::string& id);

protected:
    /*!
     *  Identify storage mode of state data (combination of properties defining state);
     *  valid modes include Phase::nativeState ("native") or other property combinations
     *  defined by Phase::fullStates (three-letter acronyms, for example "TDY", "TPX").
     */
    std::string detectMode(const std::set<std::string>& names, bool native=true);

    //! Retrieve set containing list of properties defining state
    std::set<std::string> stateProperties(const std::string& mode, bool alias=false);

    shared_ptr<Solution> m_sol; //!< Solution object associated with state data
    size_t m_size; //!< Number of entries in SolutionArray
    size_t m_stride; //!< Stride between SolutionArray entries
    AnyMap m_meta; //!< Metadata
    size_t m_index = npos; //!< Buffered index

    vector_fp m_data; //!< Work vector holding states
    std::map<std::string, vector_fp> m_extra; //!< Auxiliary data
};

}

#endif
