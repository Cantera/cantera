//! @file SolutionArray.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOLUTIONARRAY_H
#define CT_SOLUTIONARRAY_H

#include "cantera/base/global.h"
#include "cantera/base/AnyMap.h"

#if CT_USE_HIGHFIVE_HDF
namespace HighFive
{
    class File;
}
#endif

namespace Cantera
{

class Solution;
class ThermoPhase;

/*!
 *  A container class providing a convenient interface for representing many
 *  thermodynamic states using the same Solution object. C++ SolutionArray objects are
 *  one-dimensional by design; extensions to multi-dimensional arrays need to be
 *  implemented in high-level API's.
 */
class SolutionArray
{
private:
    SolutionArray(const shared_ptr<Solution>& sol,
                  size_t size,
                  const AnyMap& meta);

public:
    virtual ~SolutionArray() {}

    static shared_ptr<SolutionArray> create(const shared_ptr<Solution>& sol,
                                            size_t size=0,
                                            const AnyMap& meta={})
    {
        return shared_ptr<SolutionArray>(
            new SolutionArray(sol, size, meta));
    }

    /*!
     *  Initialize SolutionArray with independent memory management
     *
     *  @param extra Names of auxiliary data
     */
    void initialize(const std::vector<std::string>& extra={});

    /*!
     *  Size of SolutionArray (number of entries)
     */
    int size() const {
        return m_size;
    }

    /*!
     *  SolutionArray meta data.
     */
    AnyMap& meta() {
        return m_meta;
    }

    /*!
     *  Retrieve associated ThermoPhase object
     */
    std::shared_ptr<ThermoPhase> thermo();

    /*!
     *  Check whether SolutionArray contains a component.
     */
    bool hasComponent(const std::string& name) const;

    /*!
     *  Retrieve a component of the SolutionArray by name.
     */
    vector_fp getComponent(const std::string& name) const;

    /*!
     *  Set a component of the SolutionArray by name.
     *
     *  @param name  Component name
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
     *  @param root  Root location
     *  @param id  Identifier of SolutionArray root within the container file
     *  @param desc  Description
     */
    static void writeHeader(AnyMap& root, const std::string& id,
                            const std::string& desc);
#if CT_USE_HIGHFIVE_HDF
    static void writeHeader(HighFive::File& file, const std::string& id,
                            const std::string& desc);
#endif

    /*!
     *  Write SolutionArray data to container file.
     *
     *  @param root  Root location
     *  @param id  Identifier of SolutionArray within the container file
     */
    void writeEntry(AnyMap& root, const std::string& id);
#if CT_USE_HIGHFIVE_HDF
    void writeEntry(HighFive::File& file, const std::string& id);
#endif

    /*!
     *  Save the current SolutionArray to a container file.
     *
     *  @param fname  Name of output container file
     *  @param id  Identifier of SolutionArray within the container file
     */
    void save(const std::string& fname, const std::string& id, const std::string& desc);

    /*!
     *  Read header data from container file.
     *
     *  @param fname  Name of container file
     *  @param id  Identifier of SolutionArray root within the container file
     */
    static AnyMap readHeader(const std::string& fname, const std::string& id);
    static AnyMap readHeader(const AnyMap& root, const std::string& id);
#if CT_USE_HIGHFIVE_HDF
    static AnyMap readHeader(const HighFive::File& file, const std::string& id);
#endif

    /*!
     *  Restore SolutionArray from a container file.
     *
     *  @param fname  Name of container file
     *  @param id  Identifier of SolutionArray within the container file
     */
    void restore(const std::string& fname, const std::string& id);
    void restore(const AnyMap& root, const std::string& id);
#if CT_USE_HIGHFIVE_HDF
    void restore(const HighFive::File& file, const std::string& id);
#endif

protected:
    shared_ptr<Solution> m_sol; //!< Solution object associated with state data
    size_t m_size; //!< Number of entries in SolutionArray
    size_t m_stride; //!< Stride between SolutionArray entries
    AnyMap m_meta; //!< Metadata
    size_t m_index = npos; //!< Buffered index

    shared_ptr<vector_fp> m_work; //!< Work vector holding states
    double* m_data; //!< Memory location holding state information
    std::map<std::string, shared_ptr<vector_fp>> m_other; //!< Auxiliary data
};

}

#endif
