//! @file hdfUtils.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_HDF_UTILS_H
#define CT_HDF_UTILS_H

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5DataType.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

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

namespace Cantera
{

h5::Group locateH5Group(const h5::File& file, const std::string& id)
{
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    std::string grp = tokens[0];
    if (!file.exist(grp) || file.getObjectType(grp) != h5::ObjectType::Group) {
        throw CanteraError("locateH5Group", "No group with id '{}' found", grp);
    }

    std::string path = grp;
    h5::Group sub = file.getGroup(grp);
    tokens.erase(tokens.begin());
    for (auto& grp : tokens) {
        path += "/" + grp;
        if (!sub.exist(grp) || sub.getObjectType(grp) != h5::ObjectType::Group) {
            throw CanteraError("locateH5Group", "No group with id '{}' found", path);
        }
        sub = sub.getGroup(grp);
    }
    return sub;
}

h5::Group openH5Group(h5::File& file, const std::string& id)
{
    if (!file.exist(id)) {
        return file.createGroup(id);
    }
    if (file.getObjectType(id) != h5::ObjectType::Group) {
        throw CanteraError("openH5Group", "Invalid object with id '{}' exists", id);
    }
    return file.getGroup(id);
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
            double value;
            attr.read(value);
            out[name] = value;
        } else if (dclass == h5::DataTypeClass::Integer) {
            int value;
            attr.read(value);
            out[name] = value;
        } else if (dclass == h5::DataTypeClass::String) {
            std::string value;
            attr.read(value);
            out[name] = value;
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

void writeH5Attributes(h5::Group& sub, const AnyMap& meta)
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
        } else if (item.second.is<std::vector<bool>>()) {
            auto bValue = item.second.as<std::vector<bool>>();
            std::vector<H5Boolean> value;
            for (auto b : bValue) {
                value.push_back(b ? H5Boolean::TRUE : H5Boolean::FALSE);
            }
            h5::Attribute attr = sub.createAttribute<H5Boolean>(
                item.first, h5::DataSpace::From(value));
            attr.write(value);
        } else if (item.second.is<AnyMap>()) {
            // step into recursion
            auto value = item.second.as<AnyMap>();
            auto grp = sub.createGroup(item.first);
            writeH5Attributes(grp, value);
        } else {
            throw NotImplementedError("writeH5Attributes",
                "Unable to write attribute '{}' with type '{}'",
                item.first, item.second.type_str());
        }
    }
}

void writeH5FloatVector(h5::Group& sub, std::string id, vector_fp data)
{
    std::vector<size_t> dims{data.size()};
    h5::DataSet dataset = sub.createDataSet<double>(id, h5::DataSpace(dims));
    dataset.write(data);
}

vector_fp readH5FloatVector(h5::Group& sub, std::string id, size_t size)
{
    h5::DataSet dataset = sub.getDataSet(id);
    if (dataset.getDataType().getClass() != h5::DataTypeClass::Float) {
        throw CanteraError("readH5FloatVector",
            "Type of DataSet '{}' is inconsistent; expected HDF float.", id);
    }
    if (dataset.getElementCount() != size) {
        throw CanteraError("readH5FloatVector",
            "Size of DataSet '{}' is inconsistent; expected {} elements but "
            "received {} elements.", id, size, dataset.getElementCount());
    }
    vector_fp out;
    dataset.read(out);
    return out;
}

void writeH5FloatMatrix(h5::Group& sub, std::string id, std::vector<vector_fp> data)
{
    std::vector<size_t> dims{data.size()};
    dims.push_back(data.size() ? data[0].size() : 0);
    h5::DataSet dataset = sub.createDataSet<double>(id, h5::DataSpace(dims));
    dataset.write(data);
}

std::vector<vector_fp> readH5FloatMatrix(h5::Group& sub, std::string id,
                                         size_t rows, size_t cols)
{
    h5::DataSet dataset = sub.getDataSet(id);
    if (dataset.getDataType().getClass() != h5::DataTypeClass::Float) {
        throw CanteraError("readH5FloatMatrix",
            "Type of DataSet '{}' is inconsistent; expected HDF float.", id);
    }
    h5::DataSpace space = dataset.getSpace();
    if (space.getNumberDimensions() != 2) {
        throw CanteraError("readH5FloatMatrix",
            "Shape of DataSet '{}' is inconsistent; expected two dimensions.", id);
    }
    const auto& shape = space.getDimensions();
    if (shape[0] != rows) {
        throw CanteraError("readH5FloatMatrix",
            "Shape of DataSet '{}' is inconsistent; expected {} rows.", id, rows);
    }
    if (shape[1] != cols) {
        throw CanteraError("readH5FloatMatrix",
            "Shape of DataSet '{}' is inconsistent; expected {} columns.", id, cols);
    }
    std::vector<vector_fp> out;
    dataset.read(out);
    return out;
}

}

#endif
