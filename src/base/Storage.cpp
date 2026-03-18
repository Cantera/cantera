/**
 * @file Storage.cpp
 *      Definition file for class Storage.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/AnyMap.h"
#include "cantera/base/global.h"
#include "cantera/base/Storage.h"

#if CT_USE_HDF5

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5DataType.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

namespace h5 = HighFive;
typedef h5::details::Boolean H5Boolean;

#endif

namespace Cantera
{

#if CT_USE_HDF5

Storage::Storage(string fname, bool write) : m_write(write)
{
    if (m_write) {
        m_file = make_unique<h5::File>(fname, h5::File::OpenOrCreate);
    } else {
        auto fullName = findInputFile(fname);
        m_file = make_unique<h5::File>(fullName, h5::File::ReadOnly);
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

bool Storage::hasGroup(const string& id) const
{
    if (!m_file->exist(id)) {
        return false;
    }
    if (m_file->getObjectType(id) != h5::ObjectType::Group) {
        return false;
    }
    return true;
}

bool Storage::checkGroupRead(const string& id) const
{
    vector<string> tokens = tokenizePath(id);
    string grp = tokens[0];
    if (!hasGroup(grp)) {
        throw CanteraError("Storage::checkGroupRead",
             "No group with id '{}' found at root.", grp);
    }

    string path = grp;
    h5::Group sub = m_file->getGroup(grp);
    tokens.erase(tokens.begin());
    for (auto& grp : tokens) {
        if (!hasGroup(path + "/" + grp)) {
            throw CanteraError("Storage::checkGroupRead",
                "No group with id '{}' found at '{}'.", grp, path);
        }
        path += "/" + grp;
        sub = sub.getGroup(grp);
    }
    return true;
}

bool Storage::checkGroupWrite(const string& id, bool permissive)
{
    if (!m_write) {
        throw CanteraError("Storage::checkGroupWrite",
            "Cannot write to file opened in read mode.");
    }
    if (id == "") {
        throw CanteraError("Storage::checkGroupWrite",
            "Cannot write to empty group id '' (root location).");
    }
    if (!m_file->exist(id)) {
        if (!permissive) {
            throw CanteraError("Storage::checkGroupWrite",
                "Specified group with id '{}' does not exist.", id);
        }
        m_file->createGroup(id);
        return false;
    }
    if (m_file->getObjectType(id) != h5::ObjectType::Group) {
        throw CanteraError("Storage::checkGroupWrite",
            "Unable to write to existing object with id '{}'.", id);
    }
    return true;
}

bool Storage::checkGroup(const string& id, bool permissive)
{
    try {
        if (m_write) {
            return checkGroupWrite(id, permissive);
        }
        return checkGroupRead(id);
    } catch (const CanteraError& err) {
        if (permissive) {
            return false;
        }
        throw CanteraError("Storage::checkGroup", err.getMessage());
    } catch (const std::exception& err) {
        if (permissive) {
            return false;
        }
        // convert HighFive exception
        throw CanteraError("Storage::checkGroup",
            "Encountered exception for group '{}':\n{}", id, err.what());
    }
}

void Storage::deleteGroup(const string& id)
{
    try {
        m_file->unlink(id);
    } catch (const std::exception& err) {
        // convert HighFive exception
        throw CanteraError("Storage::deleteGroup",
            "Encountered exception while deleting group '{}':\n{}", id, err.what());
    }
}

pair<size_t, set<string>> Storage::contents(const string& id) const
{
    try {
        checkGroupRead(id);
    } catch (const CanteraError& err) {
        throw CanteraError("Storage::contents",
            "Caught exception for group '{}':\n", id, err.getMessage());
    }
    h5::Group sub = m_file->getGroup(id);
    set<string> names;
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
                vector<double> values;
                attr.read(values);
                out[name] = values;
            } else {
                double value;
                attr.read(value);
                out[name] = value;
            }
        } else if (dclass == h5::DataTypeClass::Integer) {
            if (attr.getSpace().getElementCount() > 1) {
                vector<long int> values;
                attr.read(values);
                out[name] = values;
            } else {
                long int value;
                attr.read(value);
                out[name] = value;
            }
        } else if (dclass == h5::DataTypeClass::String) {
            if (attr.getSpace().getElementCount() > 1) {
                vector<string> values;
                attr.read(values);
                out[name] = values;
            } else {
                string value;
                attr.read(value);
                out[name] = value;
            }
        } else if (dclass == h5::DataTypeClass::Enum) {
            // only booleans are supported
            if (attr.getSpace().getElementCount() > 1) {
                vector<H5Boolean> values;
                attr.read(values);
                vector<bool> bValues;
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

bool Storage::hasAttribute(const string& id, const string& attr) const
{
    try {
        checkGroupRead(id);
    } catch (const CanteraError& err) {
        throw CanteraError("Storage::hasAttribute",
            "Caught exception for group '{}':\n", id, err.getMessage());
    }
    h5::Group sub = m_file->getGroup(id);
    auto names = sub.listAttributeNames();
    return std::find(names.begin(), names.end(), attr) != names.end();
}

AnyMap Storage::readAttributes(const string& id, bool recursive) const
{
    try {
        checkGroupRead(id);
        h5::Group sub = m_file->getGroup(id);
        return readH5Attributes(sub, recursive);
    } catch (const Cantera::NotImplementedError& err) {
        throw NotImplementedError("Storage::readAttribute",
            "{} in group '{}'.", err.getMessage(), id);
    } catch (const CanteraError& err) {
        throw CanteraError("Storage::readAttribute",
            "Caught exception for group '{}':\n", id, err.getMessage());
    }
}

void writeH5Attributes(h5::Group sub, const AnyMap& meta)
{
    for (auto& [name, item] : meta) {
        if (sub.hasAttribute(name)) {
            throw NotImplementedError("writeH5Attributes",
                "Unable to overwrite existing Attribute '{}'", name);
        }
        if (item.is<long int>()) {
            int value = item.asInt();
            h5::Attribute attr = sub.createAttribute<long int>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<double>()) {
            double value = item.asDouble();
            h5::Attribute attr = sub.createAttribute<double>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<string>()) {
            string value = item.asString();
            h5::Attribute attr = sub.createAttribute<string>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<bool>()) {
            bool bValue = item.asBool();
            H5Boolean value = bValue ?
                H5Boolean::HighFiveTrue : H5Boolean::HighFiveFalse;
            h5::Attribute attr = sub.createAttribute<H5Boolean>(
                name, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.is<vector<long int>>()) {
            auto values = item.as<vector<long int>>();
            h5::Attribute attr = sub.createAttribute<long int>(
                name, h5::DataSpace::From(values));
            attr.write(values);
        } else if (item.is<vector<double>>()) {
            auto values = item.as<vector<double>>();
            h5::Attribute attr = sub.createAttribute<double>(
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
                values.push_back(b ?
                    H5Boolean::HighFiveTrue : H5Boolean::HighFiveFalse);
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

void Storage::writeAttributes(const string& id, const AnyMap& meta)
{
    try {
        checkGroupWrite(id, false);
        h5::Group sub = m_file->getGroup(id);
        writeH5Attributes(sub, meta);
    } catch (const Cantera::NotImplementedError& err) {
        throw NotImplementedError("Storage::writeAttribute",
            "{} in group '{}'.", err.getMessage(), id);
    } catch (const CanteraError& err) {
        // rethrow with public method attribution
        throw CanteraError("Storage::writeAttributes", "{}", err.getMessage());
    } catch (const std::exception& err) {
        // convert HighFive exception
        throw CanteraError("Storage::writeAttributes",
            "Encountered exception for group '{}':\n{}", id, err.what());
    }
}

AnyValue Storage::readData(const string& id,
                           const string& name, size_t rows, size_t cols) const
{
    try {
        checkGroupRead(id);
    } catch (const CanteraError& err) {
        throw CanteraError("Storage::readData",
            "Caught exception for group '{}':\n", id, err.getMessage());
    }
    h5::Group sub = m_file->getGroup(id);
    if (!sub.exist(name)) {
        throw CanteraError("Storage::readData",
            "DataSet '{}' not found in group '{}'.", name, id);
    }
    h5::DataSet dataset = sub.getDataSet(name);
    h5::DataSpace space = dataset.getSpace();
    size_t ndim = space.getNumberDimensions();
    if (cols == 0 && ndim != 1) {
        throw CanteraError("Storage::readData",
            "Shape of DataSet '{}' is inconsistent; expected one dimensions but "
            "received {}.", name, ndim);
    } else if (cols != 0 && cols != npos && ndim != 2) {
        throw CanteraError("Storage::readData",
            "Shape of DataSet '{}' is inconsistent; expected two dimensions but "
            "received {}.", name, ndim);
    }
    if (ndim == 0 || ndim > 2) {
        throw NotImplementedError("Storage::readData",
            "Cannot process DataSet '{}' as data has {} dimensions.", name, ndim);
    }
    const auto& shape = space.getDimensions();
    if (shape[0] != rows) {
        throw CanteraError("Storage::readData",
            "Shape of DataSet '{}' is inconsistent; expected {} rows "
            "but received {}.", name, rows, shape[0]);
    }
    if (cols != 0 && cols != npos && shape[1] != cols) {
        throw CanteraError("Storage::readData",
            "Shape of DataSet '{}' is inconsistent; expected {} columns "
            "but received {}.", name, cols, shape[1]);
    }
    AnyValue out;
    const auto datatype = dataset.getDataType().getClass();
    if (datatype == h5::DataTypeClass::Float) {
        try {
            if (ndim == 1) {
                vector<double> data;
                dataset.read(data);
                out = data;
            } else { // ndim == 2
                vector<vector<double>> data;
                dataset.read(data);
                out = data;
            }
        } catch (const std::exception& err) {
            throw NotImplementedError("Storage::readData",
                "Encountered HighFive exception for DataSet '{}' in group '{}':\n{}",
                name, id, err.what());
        }
    } else if (datatype == h5::DataTypeClass::Integer) {
        try {
            if (ndim == 1) {
                vector<long int> data;
                dataset.read(data);
                out = data;
            } else { // ndim == 2
                vector<vector<long int>> data;
                dataset.read(data);
                out = data;
            }
        } catch (const std::exception& err) {
            throw NotImplementedError("Storage::readData",
                "Encountered HighFive exception for DataSet '{}' in group '{}':\n{}",
                name, id, err.what());
        }
    } else if (datatype == h5::DataTypeClass::String) {
        try {
            if (ndim == 1) {
                vector<string> data;
                dataset.read(data);
                out = data;
            } else { // ndim == 2
                vector<vector<string>> data;
                dataset.read(data);
                out = data;
            }
        } catch (const std::exception& err) {
            throw NotImplementedError("Storage::readData",
                "Encountered HighFive exception for DataSet '{}' in group '{}':\n{}",
                name, id, err.what());
        }
    } else {
        throw NotImplementedError("Storage::readData",
            "DataSet '{}' is not readable.", name);
    }
    return out;
}

void Storage::writeData(const string& id, const string& name, const AnyValue& data)
{
    try {
        checkGroupWrite(id, false);
    } catch (const CanteraError& err) {
        // rethrow with public method attribution
        throw CanteraError("Storage::writeData", "{}", err.getMessage());
    } catch (const std::exception& err) {
        // convert HighFive exception
        throw CanteraError("Storage::writeData",
            "Encountered exception for group '{}':\n{}", id, err.what());
    }
    h5::Group sub = m_file->getGroup(id);
    if (sub.exist(name)) {
        throw NotImplementedError("Storage::writeData",
            "Unable to overwrite existing DataSet '{}' in group '{}'.", name, id);
    }
    size_t size = data.vectorSize();
    auto [rows, cols] = data.matrixShape();
    if (size == npos && rows == npos) {
        throw CanteraError("Storage::writeData",
            "Cannot write DataSet '{}' in group '{}' as input data with type\n"
            "'{}'\nis neither a vector nor a matrix.", name, id, data.type_str());
    }
    vector<size_t> dims{data.vectorSize()};
    if (data.isVector<long int>()) {
        h5::DataSet dataset = sub.createDataSet<long int>(name, h5::DataSpace(dims));
        dataset.write(data.asVector<long int>());
        return;
    }
    if (data.isVector<double>()) {
        h5::DataSet dataset = sub.createDataSet<double>(name, h5::DataSpace(dims));
        dataset.write(data.asVector<double>());
        return;
    }
    if (data.isVector<string>()) {
        h5::DataSet dataset = sub.createDataSet<string>(name, h5::DataSpace(dims));
        dataset.write(data.asVector<string>());
        return;
    }
    if (cols != npos) {
        dims.clear();
        dims.push_back(rows);
        dims.push_back(cols);
    } else {
        throw NotImplementedError("Storage::writeData",
            "Cannot write DataSet '{}' in group '{}' as input data with type\n"
            "'{}'\nis not supported.", name, id, data.type_str());
    }
    if (m_compressionLevel) {
        // Set chunk size to single chunk and apply compression level; for caveats, see
        // https://stackoverflow.com/questions/32994766/compressed-files-bigger-in-h5py
        h5::DataSpace space(dims, dims); //{h5::DataSpace::UNLIMITED, dims[1]});
        h5::DataSetCreateProps props;
        props.add(h5::Chunking(vector<hsize_t>{dims[0], dims[1]}));
        props.add(h5::Deflate(m_compressionLevel));
        if (data.isVector<vector<long int>>()) {
            h5::DataSet dataset = sub.createDataSet<long int>(name, space, props);
            dataset.write(data.asVector<vector<long int>>());
        } else if (data.isVector<vector<double>>()) {
            h5::DataSet dataset = sub.createDataSet<double>(name, space, props);
            dataset.write(data.asVector<vector<double>>());
        } else if (data.isVector<vector<string>>()) {
            h5::DataSet dataset = sub.createDataSet<string>(name, space, props);
            dataset.write(data.asVector<vector<string>>());
        } else {
            throw NotImplementedError("Storage::writeData",
                "Cannot write DataSet '{}' in group '{}' as input data with type\n"
                "'{}'\nis not supported.", name, id, data.type_str());
        }
    } else {
        h5::DataSpace space(dims);
        if (data.isVector<vector<long int>>()) {
            h5::DataSet dataset = sub.createDataSet<long int>(name, space);
            dataset.write(data.asVector<vector<long int>>());
        } else if (data.isVector<vector<double>>()) {
            h5::DataSet dataset = sub.createDataSet<double>(name, space);
            dataset.write(data.asVector<vector<double>>());
        } else if (data.isVector<vector<string>>()) {
            h5::DataSet dataset = sub.createDataSet<string>(name, space);
            dataset.write(data.asVector<vector<string>>());
        } else {
            throw NotImplementedError("Storage::writeData",
                "Cannot write DataSet '{}' in group '{}' as input data with type\n"
                "'{}'\nis not supported.", name, id, data.type_str());
        }
    }
}

#else

Storage::Storage(string fname, bool write)
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

bool Storage::hasGroup(const string& id) const
{
    throw CanteraError("Storage::hasGroup",
                       "Saving to HDF requires HighFive installation.");
}

bool Storage::checkGroup(const string& id, bool permissive)
{
    throw CanteraError("Storage::checkGroup",
                       "Saving to HDF requires HighFive installation.");
}

void Storage::deleteGroup(const string& id)
{
    throw CanteraError("Storage::deleteGroup",
                       "Saving to HDF requires HighFive installation.");
}

pair<size_t, set<string>> Storage::contents(const string& id) const
{
    throw CanteraError("Storage::contents",
                       "Saving to HDF requires HighFive installation.");
}

bool Storage::hasAttribute(const string& id, const string& attr) const
{
    throw CanteraError("Storage::hasAttribute",
                       "Saving to HDF requires HighFive installation.");
}

AnyMap Storage::readAttributes(const string& id, bool recursive) const
{
    throw CanteraError("Storage::readAttributes",
                       "Saving to HDF requires HighFive installation.");
}

void Storage::writeAttributes(const string& id, const AnyMap& meta)
{
    throw CanteraError("Storage::writeAttributes",
                       "Saving to HDF requires HighFive installation.");
}

AnyValue Storage::readData(const string& id,
                           const string& name, size_t rows, size_t cols) const
{
    throw CanteraError("Storage::readData",
                       "Saving to HDF requires HighFive installation.");
}

void Storage::writeData(const string& id,
                        const string& name, const AnyValue& data)
{
    throw CanteraError("Storage::writeData",
                       "Saving to HDF requires HighFive installation.");
}

#endif

}
