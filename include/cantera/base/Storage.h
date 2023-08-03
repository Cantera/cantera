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

/**
 *  A wrapper class handling storage to HDF. Acts as a thin wrapper for HighFive.
 *  The class implements methods that are intended to be called from SolutionArray.
 *
 *  @since New in %Cantera 3.0.
 *  @warning This class is an experimental part of the %Cantera API and may be
 *      changed or removed without notice.
 */
class Storage
{
public:
    Storage(string fname, bool write);

    ~Storage();

    //! Set compression level (0..9)
    //!
    //! Compression is only applied to matrix-type data; note that compression may
    //! increase file size for small data sets (compression requires setting of chunk
    //! sizes, which involves considerable overhead for metadata).
    void setCompressionLevel(int level);

    //! Check whether location `id` represents a group
    bool hasGroup(const string& id) const;

    //! Check whether path location exists.
    //! If the location does not exist, an exception is thrown unless the *permissive*
    //! flag is set; in this case, the method attempts to create a new location if the
    //! file is accessed in write mode.
    //! @param id  storage location within file
    //! @param permissive  if true, do not raise exceptions (default=false)
    //! @returns  boolean indicating whether id is pre-existing
    bool checkGroup(const string& id, bool permissive=false);

    //! Delete group
    //! @param id  storage location within file
    void deleteGroup(const string& id);

    //! Retrieve contents of file from a specified location
    //! @param id  storage location within file
    //! @returns  pair containing size and list of entry names of stored data set
    pair<size_t, set<string>> contents(const string& id) const;

    //! Read attributes from a specified location
    //! @param id  storage location within file
    //! @param attr  name of attribute to be checked
    bool hasAttribute(const string& id, const string& attr) const;

    //! Read attributes from a specified location
    //! @param id  storage location within file
    //! @param recursive  boolean indicating whether subgroups should be included
    //! @returns  AnyMap containing attributes
    AnyMap readAttributes(const string& id, bool recursive) const;

    //! Write attributes to a specified location
    //! @param id  storage location within file
    //! @param meta  AnyMap containing attributes
    void writeAttributes(const string& id, const AnyMap& meta);

    //! Read dataset from a specified location
    //! @param id  storage location within file
    //! @param name  name of vector/matrix entry
    //! @param rows  number of vector length or matrix rows
    //! @param cols  number of matrix columns, if applicable; if 0, a vector is
    //!     expected, if npos, the size is detected automatically; otherwise, an exact
    //!     number of columns needs to be matched.
    //! @returns  matrix or vector containing data; implemented for types
    //!     `vector<double>`, `vector<long int>`, `vector<string>`,
    //!     `vector<vector<double>>`, `vector<vector<long int>>` and
    //!     `vector<vector<string>>`
    AnyValue readData(const string& id,
                      const string& name, size_t rows, size_t cols=npos) const;

    //! Write dataset to a specified location
    //! @param id  storage location within file
    //! @param name  name of matrix entry
    //! @param data  vector or matrix containing data; implemented for types
    //!     `vector<double>`, `vector<long int>`, `vector<string>`
    //!     `vector<vector<double>>`, `vector<vector<long int>>` and
    //!     `vector<vector<string>>`
    void writeData(const string& id, const string& name, const AnyValue& data);

private:
#if CT_USE_HDF5
    //! ensure that HDF group is readable
    bool checkGroupRead(const string& id) const;

    //! ensure that HDF group is writeable
    bool checkGroupWrite(const string& id, bool permissive);

    unique_ptr<HighFive::File> m_file; //!< HDF container file
    bool m_write; //!< HDF access mode
    int m_compressionLevel=0; //!< HDF compression level
#endif
};

}

#endif
