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

/**
 * A container class providing a convenient interface for representing many
 * thermodynamic states using the same Solution object. C++ SolutionArray objects are
 * one-dimensional by design; extensions to multi-dimensional arrays need to be
 * implemented in high-level API's.
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

    /**
     * Initialize SolutionArray with independent memory management
     *
     * @param extra Names of auxiliary data
     */
    void initialize(const std::vector<std::string>& extra={});

    /**
     * Initialize SolutionArray object with mapped memory
     *
     * @param data Pointer to mapped memory address
     * @param size Number of entries in SolutionArray
     * @param stride An integer indicating the stride between entries
     * @param offsets A vector of pairs containing offsets within the mapped memory
     */
    void initialize(double* data,
                    size_t size,
                    size_t stride,
                    const std::vector<std::pair<std::string, size_t>>& offsets);

    /**
     * Size of SolutionArray (number of entries)
     */
    int size() const {
        return m_size;
    }

    /**
     * Save the current SolutionArray to a container file.
     *
     * @param fname Name of output container file
     * @param id Identifier of SolutionArray within the container file
     */
    void save(const std::string& fname, const std::string& id);

    /**
     * Restore SolutionArray from a container file.
     *
     * @param fname Name of container file
     * @param id Identifier of SolutionArray within the container file
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
    bool m_managed = false; //!< Flag indicating whether memory is externally managed

    shared_ptr<vector_fp> m_work; //!< Work vector holding states (if not managed)
    double* m_data; //!< Memory location holding state information (may be augmented)
    std::map<std::string, shared_ptr<vector_fp>> m_other; //!< Auxiliary data
    std::map<std::string, size_t> m_offsets; //!< Map of offsets in state vector
    std::map<std::string, size_t> m_extra; //!< Map of offsets in auxiliary data
};


// /**
//  * Create a SolutionArray object with independent memory management
//  *
//  * @param sol The Solution object associated with state information
//  * @param max_size Expected maximum number of entries in SolutionArray
//  * @param extra A vector of additional entries
//  * @param meta Metadata
//  * @return shared_ptr<SolutionArray>
//  */
// shared_ptr<SolutionArray> newSolutionArray(
//     const shared_ptr<Solution>& sol,
//     size_t max_size,
//     const std::vector<std::string>& extra={},
//     const AnyMap& meta={})
// {
//     shared_ptr<SolutionArray> arr = SolutionArray::create(sol, max_size, meta);
//     arr->initialize(extra);
//     return arr;
// }

}

#endif
