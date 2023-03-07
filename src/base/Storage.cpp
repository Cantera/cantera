/**
 * @file Storage.cpp
 *      Definition file for class Storage.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/AnyMap.h"
#include "cantera/base/Storage.h"

#if CT_USE_HDF5

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

#if CT_USE_HDF5

Storage::Storage(std::string fname, bool write) : m_write(write)
{
    if (m_write) {
        m_file = make_unique<h5::File>(fname, h5::File::OpenOrCreate);
    } else {
        m_file = make_unique<h5::File>(fname, h5::File::ReadOnly);
    }
}

Storage::~Storage()
{
    m_file->flush();
}

void Storage::setCompressionLevel(int level)
{
    if (level < 0 || level > 9) {
        throw CanteraError("Storage::setCompressionLevel",
                           "Invalid compression level '{}' (needs to be 0..9).", level);
    }
    m_compressionLevel = level;
}

bool Storage::checkGroupRead(const std::string& id) const
{
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    std::string grp = tokens[0];
    if (!m_file->exist(grp) || m_file->getObjectType(grp) != h5::ObjectType::Group) {
        throw CanteraError("Storage::checkGroupRead",
                           "No group with id '{}' found", grp);
    }

    std::string path = grp;
    h5::Group sub = m_file->getGroup(grp);
    tokens.erase(tokens.begin());
    for (auto& grp : tokens) {
        path += "/" + grp;
        if (!sub.exist(grp) || sub.getObjectType(grp) != h5::ObjectType::Group) {
            throw CanteraError("Storage::checkGroupRead",
                               "No group with id '{}' found", path);
        }
        sub = sub.getGroup(grp);
    }
    return true;
}

bool Storage::checkGroupWrite(const std::string& id)
{
    if (!m_file->exist(id)) {
        m_file->createGroup(id);
        return true;
    }
    if (m_file->getObjectType(id) != h5::ObjectType::Group) {
        throw CanteraError("Storage::checkGroupWrite",
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
    h5::Group sub = m_file->getGroup(id);
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
            "Encountered invalid data with {} dimensions.", nDims);
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
    h5::Group sub = m_file->getGroup(id);
    try {
        return readH5Attributes(sub, recursive);
    } catch (const Cantera::NotImplementedError& err) {
        throw NotImplementedError("Storage::readAttribute",
            "{} in group '{}'.", err.getMessage(), id);
    }
}

void writeH5Attributes(h5::Group sub, const AnyMap& meta)
{
    for (auto& [name, item] : meta) {
        if (sub.hasAttribute(name)) {
            throw NotImplementedError("writeH5Attributes",
                "Unable to overwrite existing Attribute '{}'", name);
        }
        if (item.is<double>()) {
            double value = item.asDouble();
            h5::Attribute attr = sub.createAttribute<double>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<int>() || item.is<long int>()) {
            int value = item.asInt();
            h5::Attribute attr = sub.createAttribute<int>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<string>()) {
            string value = item.asString();
            h5::Attribute attr = sub.createAttribute<string>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<bool>()) {
            bool bValue = item.asBool();
            H5Boolean value = bValue ? H5Boolean::TRUE : H5Boolean::FALSE;
            h5::Attribute attr = sub.createAttribute<H5Boolean>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<vector<double>>()) {
            auto values = item.as<vector<double>>();
            h5::Attribute attr = sub.createAttribute<double>(
                name, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.is<vector<int>>()) {
            auto values = item.as<vector<int>>();
            h5::Attribute attr = sub.createAttribute<int>(
                name, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.is<vector<string>>()) {
            auto values = item.as<vector<string>>();
            h5::Attribute attr = sub.createAttribute<string>(
                name, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.is<vector<bool>>()) {
            auto bValue = item.as<vector<bool>>();
            vector<H5Boolean> values;
            for (auto b : bValue) {
                values.push_back(b ? H5Boolean::TRUE : H5Boolean::FALSE);
            }
            h5::Attribute attr = sub.createAttribute<H5Boolean>(
                name, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.is<AnyMap>()) {
            // step into recursion
            auto value = item.as<AnyMap>();
            auto grp = sub.createGroup(name);
            writeH5Attributes(grp, value);
        } else {
            throw NotImplementedError("writeH5Attributes",
                "Unable to write attribute '{}' with type '{}'",
                name, item.type_str());
        }
    }
}

void Storage::writeAttributes(const std::string& id, const AnyMap& meta)
{
    h5::Group sub = m_file->getGroup(id);
    try {
        writeH5Attributes(sub, meta);
    } catch (const Cantera::NotImplementedError& err) {
        throw NotImplementedError("Storage::writeAttribute",
            "{} in group '{}'.", err.getMessage(), id);
    }
}

vector_fp Storage::readVector(const std::string& id,
                              const std::string& name, size_t size) const
{
    h5::Group sub = m_file->getGroup(id);
    if (!sub.exist(name)) {
        throw CanteraError("Storage::readVector",
            "DataSet '{}' not found in group '{}'.", name, id);
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
    h5::Group sub = m_file->getGroup(id);
    if (sub.exist(name)) {
        throw NotImplementedError("Storage::writeVector",
            "Unable to overwrite existing DataSet '{}' in group '{}'.", name, id);
    }
    std::vector<size_t> dims{data.size()};
    h5::DataSet dataset = sub.createDataSet<double>(name, h5::DataSpace(dims));
    dataset.write(data);
}

std::vector<vector_fp> Storage::readMatrix(const std::string& id,
                                           const std::string& name,
                                           size_t rows, size_t cols) const
{
    h5::Group sub = m_file->getGroup(id);
    if (!sub.exist(name)) {
        throw CanteraError("Storage::readMatrix",
            "DataSet '{}' not found in group '{}'.", name, id);
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
    h5::Group sub = m_file->getGroup(id);
    if (sub.exist(name)) {
        throw NotImplementedError("Storage::writeMatrix",
            "Unable to overwrite existing DataSet '{}' in group '{}'.", name, id);
    }
    std::vector<size_t> dims{data.size()};
    dims.push_back(data.size() ? data[0].size() : 0);
    if (m_compressionLevel) {
        // Set chunk size to single chunk and apply compression level; for caveats, see
        // https://stackoverflow.com/questions/32994766/compressed-files-bigger-in-h5py
        h5::DataSpace space(dims, dims); //{h5::DataSpace::UNLIMITED, dims[1]});
        h5::DataSetCreateProps props;
        props.add(h5::Chunking(std::vector<hsize_t>{dims[0], dims[1]}));
        props.add(h5::Deflate(m_compressionLevel));
        h5::DataSet dataset = sub.createDataSet<double>(name, space, props);
        dataset.write(data);
    } else {
        h5::DataSpace space(dims);
        h5::DataSet dataset = sub.createDataSet<double>(name, space);
        dataset.write(data);
    }
}

#else

Storage::Storage(std::string fname, bool write)
{
    throw CanteraError("Storage::Storage",
                       "Saving to HDF requires HighFive installation.");
}

Storage::~Storage()
{
}

void Storage::setCompressionLevel(int level)
{
    throw CanteraError("Storage::setCompressionLevel",
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
