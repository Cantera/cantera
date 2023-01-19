//! @file Storage.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STORAGE_H
#define CT_STORAGE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/stringUtils.h"

#if CT_USE_HDF5

#ifdef _WIN32
  // see https://github.com/microsoft/vcpkg/issues/24293
  #define H5_BUILT_AS_DYNAMIC_LIB
#else
  #define H5_BUILT_AS_STATIC_LIB
#endif

namespace HighFive {
    class File;
}

#endif

namespace Cantera
{

/*!
 *  A wrapper class handling storage to HDF; acts as a thin wrapper for HighFive
 *
 *  @since  New in Cantera 3.0.
 *  @warning This class is an experimental part of the %Cantera API and may be
 *      changed or removed without notice.
 */
class Storage
{
public:
    Storage(std::string fname, bool write);

    ~Storage();

    //! Set compression level (0..9)
    /*!
     *  Compression is only applied to species data; note that compression may increase
     *  file size for small data sets (compression requires setting of chunk sizes,
     *  which involves considerable overhead for metadata).
     */
    void setCompressionLevel(int level);

    //! Check whether path location exists
    //! If the file has write access, create location if necessary
    //! @param id  storage location within file
    bool checkGroup(const std::string& id);

    //! Retrieve contents of file from a specified location
    //! @param id  storage location within file
    //! @returns  pair containing size and list of entry names of stored data set
    std::pair<size_t, std::set<std::string>> contents(const std::string& id) const;

    //! Read attributes from a specified location
    //! @param id  storage location within file
    //! @param attr  name of attribute to be checked
    bool hasAttribute(const std::string& id, const std::string& attr) const;

    //! Read attributes from a specified location
    //! @param id  storage location within file
    //! @param recursive  boolean indicating whether subgroups should be included
    //! @returns  AnyMap containing attributes
    AnyMap readAttributes(const std::string& id, bool recursive) const;

    //! Write attributes to a specified location
    //! @param id  storage location within file
    //! @param meta  AnyMap containing attributes
    void writeAttributes(const std::string& id, const AnyMap& meta);

    //! Read data vector from a specified location
    //! @param id  storage location within file
    //! @param name  name of data vector entry
    //! @param size  size of data vector entry
    //! @returns  data vector
    vector_fp readVector(const std::string& id,
                         const std::string& name, size_t size) const;

    //! Write data vector to a specified location
    //! @param id  storage location within file
    //! @param name  name of data vector entry
    //! @param data  data vector
    void writeVector(const std::string& id,
                     const std::string& name, const vector_fp& data);

    //! Read matrix from a specified location
    //! @param id  storage location within file
    //! @param name  name of matrix entry
    //! @param rows  number of matrix rows
    //! @param cols  number of matrix columns
    //! @returns  matrix containing data (vector of vectors)
    std::vector<vector_fp> readMatrix(const std::string& id,
                                      const std::string& name,
                                      size_t rows, size_t cols) const;

    //! Write matrix to a specified location
    //! @param id  storage location within file
    //! @param name  name of matrix entry
    //! @param data  matrix containing data (vector of vectors)
    void writeMatrix(const std::string& id,
                     const std::string& name, const std::vector<vector_fp>& data);

private:
#if CT_USE_HDF5
    bool checkGroupRead(const std::string& id) const;
    bool checkGroupWrite(const std::string& id);

    std::unique_ptr<HighFive::File> m_file;
    bool m_write;
    int m_compressionLevel=0;
#endif
};

}

#endif
