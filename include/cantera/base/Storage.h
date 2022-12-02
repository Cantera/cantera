//! @file Storage.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STORAGE_H
#define CT_STORAGE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/stringUtils.h"
#include <set>

#if CT_USE_HIGHFIVE_HDF
#if CT_USE_SYSTEM_HIGHFIVE
  #include <highfive/H5Attribute.hpp>
  #include <highfive/H5DataSet.hpp>
  #include <highfive/H5DataSpace.hpp>
  #include <highfive/H5DataType.hpp>
  #include <highfive/H5File.hpp>
  #include <highfive/H5Group.hpp>
#else
  #include "cantera/ext/HighFive/H5Attribute.hpp"
  #include "cantera/ext/HighFive/H5DataSet.hpp"
  #include "cantera/ext/HighFive/H5DataSpace.hpp"
  #include "cantera/ext/HighFive/H5DataType.hpp"
  #include "cantera/ext/HighFive/H5File.hpp"
  #include "cantera/ext/HighFive/H5Group.hpp"
#endif

namespace h5 = HighFive;

enum class H5Boolean {
    FALSE = 0,
    TRUE = 1,
};

h5::EnumType<H5Boolean> create_enum_boolean() {
    return {{"FALSE", H5Boolean::FALSE},
            {"TRUE", H5Boolean::TRUE}};
}

HIGHFIVE_REGISTER_TYPE(H5Boolean, create_enum_boolean)

#endif

namespace Cantera
{

/*!
 *  A wrapper class handling storage to HDF; acts as a thin wrapper for HighFive
 */
class Storage
{
public:
#if CT_USE_HIGHFIVE_HDF
    Storage(h5::File file, bool write) : m_file(file), m_write(write) {}
#else
    Storage() {
        throw CanteraError("Storage::Storage",
                           "Instantiation of Storage requires HighFive::File object.");
    }
#endif

    //! Flush file contents
    void flush();

    //! Check whether path go location exists
    //! If the file has write access, create location if necessary
    bool checkGroup(const std::string& id);

    //! Retrieve contents of file from a specified location
    std::pair<size_t, std::set<std::string>> contents(const std::string& id) const;

    //! Read attributes from a specified location
    AnyMap readAttributes(const std::string& id, bool recursive) const;

    //! Write attributes to a specified location
    void writeAttributes(const std::string& id, const AnyMap& meta);

    //! Read data vector from a specified location
    vector_fp readVector(const std::string& id,
                         const std::string& name, size_t size) const;

    //! Write data vector to a specified location
    void writeVector(const std::string& id,
                     const std::string& name, const vector_fp& data);

    //! Read matrix from a specified location
    std::vector<vector_fp> readMatrix(const std::string& id,
                                      const std::string& name,
                                      size_t rows, size_t cols) const;

    //! Write matrix to a specified location
    void writeMatrix(const std::string& id,
                     const std::string& name, const std::vector<vector_fp>& data);

private:
#if CT_USE_HIGHFIVE_HDF
    bool checkGroupRead(const std::string& id) const;
    bool checkGroupWrite(const std::string& id);

    h5::File m_file;
#endif

    bool m_write;
};

#if CT_USE_HIGHFIVE_HDF

void Storage::flush()
{
    m_file.flush();
}

bool Storage::checkGroupRead(const std::string& id) const
{
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    std::string grp = tokens[0];
    if (!m_file.exist(grp) || m_file.getObjectType(grp) != h5::ObjectType::Group) {
        throw CanteraError("Storage::checkGroup",
                           "No group with id '{}' found", grp);
    }

    std::string path = grp;
    h5::Group sub = m_file.getGroup(grp);
    tokens.erase(tokens.begin());
    for (auto& grp : tokens) {
        path += "/" + grp;
        if (!sub.exist(grp) || sub.getObjectType(grp) != h5::ObjectType::Group) {
            throw CanteraError("Storage::checkGroup",
                               "No group with id '{}' found", path);
        }
        sub = sub.getGroup(grp);
    }
    return true;
}

bool Storage::checkGroupWrite(const std::string& id)
{
    if (!m_file.exist(id)) {
        m_file.createGroup(id);
        return true;
    }
    if (m_file.getObjectType(id) != h5::ObjectType::Group) {
        throw CanteraError("Storage::checkGroup",
                           "Invalid object with id '{}' exists", id);
    }
    return true;
}

bool Storage::checkGroup(const std::string& id) {
    if (m_write) {
        return checkGroupWrite(id);
    }
    return checkGroupRead(id);
}

std::pair<size_t, std::set<std::string>> Storage::contents(const std::string& id) const
{
    h5::Group sub = m_file.getGroup(id);
    std::set<std::string> names;
    size_t nDims = npos;
    size_t nElements = 0;
    for (auto& name : sub.listObjectNames()) {
        if (sub.getObjectType(name) == h5::ObjectType::Dataset) {
            h5::DataSpace space = sub.getDataSet(name).getSpace();
            names.insert(name);
            if (space.getNumberDimensions() < nDims) {
                nDims = space.getNumberDimensions();
                nElements = space.getElementCount();
            }
        }
    }
    if (nDims != 1 && nDims != npos) {
        throw NotImplementedError("Storage::content",
            "Unable to restore data with {} dimensions.", nDims);
    }
    return std::make_pair(nElements, names);
}

AnyMap readH5Attributes(const h5::Group& sub, bool recursive)
{
    // restore meta data from attributes
    AnyMap out;
    for (auto& name : sub.listAttributeNames()) {
        h5::Attribute attr = sub.getAttribute(name);
        h5::DataType dtype = attr.getDataType();
        h5::DataTypeClass dclass = dtype.getClass();
        if (dclass == h5::DataTypeClass::Float) {
            if (attr.getSpace().getElementCount() > 1) {
                std::vector<double> values;
                attr.read(values);
                out[name] = values;
            } else {
                double value;
                attr.read(value);
                out[name] = value;
            }
        } else if (dclass == h5::DataTypeClass::Integer) {
            if (attr.getSpace().getElementCount() > 1) {
                std::vector<int> values;
                attr.read(values);
                out[name] = values;
            } else {
                int value;
                attr.read(value);
                out[name] = value;
            }
        } else if (dclass == h5::DataTypeClass::String) {
            if (attr.getSpace().getElementCount() > 1) {
                std::vector<std::string> values;
                attr.read(values);
                out[name] = values;
            } else {
                std::string value;
                attr.read(value);
                out[name] = value;
            }
        } else if (dclass == h5::DataTypeClass::Enum) {
            // only booleans are supported
            if (attr.getSpace().getElementCount() > 1) {
                std::vector<H5Boolean> values;
                attr.read(values);
                std::vector<bool> bValues;
                for (auto v : values) {
                    bValues.push_back(bool(v));
                }
                out[name] = bValues;
            } else {
                H5Boolean value;
                attr.read(value);
                out[name] = bool(value);
            }
        } else {
            throw NotImplementedError("readH5Attributes",
                "Unable to read attribute '{}' with type '{}'", name, dtype.string());
        }
    }

    if (recursive) {
        for (auto& name : sub.listObjectNames()) {
            if (sub.getObjectType(name) == h5::ObjectType::Group) {
                out[name] = readH5Attributes(sub.getGroup(name), recursive);
            }
        }
    }

    return out;
}

AnyMap Storage::readAttributes(const std::string& id, bool recursive) const
{
    h5::Group sub = m_file.getGroup(id);
    return readH5Attributes(sub, recursive);
}

void writeH5Attributes(h5::Group sub, const AnyMap& meta)
{
    for (auto& item : meta) {
        if (item.second.is<double>()) {
            double value = item.second.asDouble();
            h5::Attribute attr = sub.createAttribute<double>(
                item.first, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.second.is<int>() || item.second.is<long int>()) {
            int value = item.second.asInt();
            h5::Attribute attr = sub.createAttribute<int>(
                item.first, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.second.is<std::string>()) {
            std::string value = item.second.asString();
            h5::Attribute attr = sub.createAttribute<std::string>(
                item.first, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.second.is<bool>()) {
            bool bValue = item.second.asBool();
            H5Boolean value = bValue ? H5Boolean::TRUE : H5Boolean::FALSE;
            h5::Attribute attr = sub.createAttribute<H5Boolean>(
                item.first, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.second.is<std::vector<double>>()) {
            auto values = item.second.as<std::vector<double>>();
            h5::Attribute attr = sub.createAttribute<double>(
                item.first, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.second.is<std::vector<int>>()) {
            auto values = item.second.as<std::vector<int>>();
            h5::Attribute attr = sub.createAttribute<int>(
                item.first, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.second.is<std::vector<std::string>>()) {
            auto values = item.second.as<std::vector<std::string>>();
            h5::Attribute attr = sub.createAttribute<std::string>(
                item.first, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.second.is<std::vector<bool>>()) {
            auto bValue = item.second.as<std::vector<bool>>();
            std::vector<H5Boolean> values;
            for (auto b : bValue) {
                values.push_back(b ? H5Boolean::TRUE : H5Boolean::FALSE);
            }
            h5::Attribute attr = sub.createAttribute<H5Boolean>(
                item.first, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.second.is<AnyMap>()) {
            // step into recursion
            auto value = item.second.as<AnyMap>();
            auto grp = sub.createGroup(item.first);
            writeH5Attributes(grp, value);
        } else {
            throw NotImplementedError("Storage::writeAttributes",
                "Unable to write attribute '{}' with type '{}'",
                item.first, item.second.type_str());
        }
    }
}

void Storage::writeAttributes(const std::string& id, const AnyMap& meta)
{
    h5::Group sub = m_file.getGroup(id);
    writeH5Attributes(sub, meta);
}

vector_fp Storage::readVector(const std::string& id,
                              const std::string& name, size_t size) const
{
    h5::Group sub = m_file.getGroup(id);
    if (!sub.exist(name)) {
        throw CanteraError("Storage::readVector",
            "DataSet '{}' not found in path '{}'.", name, id);
    }
    h5::DataSet dataset = sub.getDataSet(name);
    if (dataset.getDataType().getClass() != h5::DataTypeClass::Float) {
        throw CanteraError("Storage::readVector",
            "Type of DataSet '{}' is inconsistent; expected HDF float.", name);
    }
    if (dataset.getElementCount() != size) {
        throw CanteraError("Storage::readVector",
            "Size of DataSet '{}' is inconsistent; expected {} elements but "
            "received {} elements.", name, size, dataset.getElementCount());
    }
    vector_fp out;
    dataset.read(out);
    return out;
}

void Storage::writeVector(const std::string& id,
                          const std::string& name, const vector_fp& data)
{
    h5::Group sub = m_file.getGroup(id);
    std::vector<size_t> dims{data.size()};
    h5::DataSet dataset = sub.createDataSet<double>(name, h5::DataSpace(dims));
    dataset.write(data);
}

std::vector<vector_fp> Storage::readMatrix(const std::string& id,
                                           const std::string& name,
                                           size_t rows, size_t cols) const
{
    h5::Group sub = m_file.getGroup(id);
    if (!sub.exist(name)) {
        throw CanteraError("Storage::readVector",
            "DataSet '{}' not found in path '{}'.", name, id);
    }
    h5::DataSet dataset = sub.getDataSet(name);
    if (dataset.getDataType().getClass() != h5::DataTypeClass::Float) {
        throw CanteraError("Storage::readMatrix",
            "Type of DataSet '{}' is inconsistent; expected HDF float.", name);
    }
    h5::DataSpace space = dataset.getSpace();
    if (space.getNumberDimensions() != 2) {
        throw CanteraError("Storage::readMatrix",
            "Shape of DataSet '{}' is inconsistent; expected two dimensions.", name);
    }
    const auto& shape = space.getDimensions();
    if (shape[0] != rows) {
        throw CanteraError("Storage::readMatrix",
            "Shape of DataSet '{}' is inconsistent; expected {} rows.", name, rows);
    }
    if (shape[1] != cols) {
        throw CanteraError("Storage::readMatrix",
            "Shape of DataSet '{}' is inconsistent; expected {} columns.", name, cols);
    }
    std::vector<vector_fp> out;
    dataset.read(out);
    return out;
}

void Storage::writeMatrix(const std::string& id,
                          const std::string& name, const std::vector<vector_fp>& data)
{
    h5::Group sub = m_file.getGroup(id);
    std::vector<size_t> dims{data.size()};
    dims.push_back(data.size() ? data[0].size() : 0);
    h5::DataSet dataset = sub.createDataSet<double>(name, h5::DataSpace(dims));
    dataset.write(data);
}

#else

void Storage::flush()
{
    throw CanteraError("Storage::flush",
                        "Saving to HDF requires HighFive installation.");
}

bool Storage::checkGroup(const std::string& id)
{
    throw CanteraError("Storage::checkGroup",
                        "Saving to HDF requires HighFive installation.");
}

std::pair<size_t, std::set<std::string>> Storage::contents(const std::string& id) const
{
    throw CanteraError("Storage::contents",
                        "Saving to HDF requires HighFive installation.");
}

AnyMap Storage::readAttributes(const std::string& id, bool recursive) const
{
    throw CanteraError("Storage::readAttributes",
                        "Saving to HDF requires HighFive installation.");
}

void Storage::writeAttributes(const std::string& id, const AnyMap& meta)
{
    throw CanteraError("Storage::writeAttributes",
                        "Saving to HDF requires HighFive installation.");
}

vector_fp Storage::readVector(const std::string& id,
                              const std::string& name, size_t size) const
{
    throw CanteraError("Storage::readVector",
                        "Saving to HDF requires HighFive installation.");
}

void Storage::writeVector(const std::string& id,
                          const std::string& name, const vector_fp& data)
{
    throw CanteraError("Storage::writeVector",
                        "Saving to HDF requires HighFive installation.");
}

std::vector<vector_fp> Storage::readMatrix(const std::string& id,
                                           const std::string& name,
                                           size_t rows, size_t cols) const
{
    throw CanteraError("Storage::readMatrix",
                        "Saving to HDF requires HighFive installation.");
}

void Storage::writeMatrix(const std::string& id,
                          const std::string& name, const std::vector<vector_fp>& data)
{
    throw CanteraError("Storage::writeMatrix",
                        "Saving to HDF requires HighFive installation.");
}

#endif

}

#endif
