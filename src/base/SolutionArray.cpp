/**
 * @file SolutionArray.cpp
 *      Definition file for class SolutionArray.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/SolutionArray.h"
#include "cantera/base/Solution.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include <set>

#if CT_USE_HIGHFIVE_HDF
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
#endif

namespace Cantera
{

SolutionArray::SolutionArray(
    const shared_ptr<Solution>& sol,
    size_t size,
    const AnyMap& meta)
    : m_sol(sol)
    , m_size(size)
    , m_meta(meta)
{
    if (!m_sol) {
        throw CanteraError("SolutionArray::SolutionArray",
            "Unable to create SolutionArray from invalid Solution object.");
    }
}

void SolutionArray::initialize(const std::vector<std::string>& extra)
{
    m_stride = m_sol->thermo()->stateSize();
    m_work.reset(new vector_fp(m_size * m_stride, 0.));
    m_data = m_work->data();
    for (auto& key : extra) {
        m_other.emplace(key, std::make_shared<vector_fp>(m_size));
    }
}

std::shared_ptr<ThermoPhase> SolutionArray::thermo()
{
    return m_sol->thermo();
}

bool SolutionArray::hasComponent(const std::string& name) const
{
    if (m_other.count(name)) {
        // auxiliary data
        return true;
    }
    if (m_sol->thermo()->speciesIndex(name) != npos) {
        // species
        return true;
    }
    if (name == "X" || name == "Y") {
        // reserved names
        return false;
    }
    // native state
    return (m_sol->thermo()->nativeState().count(name));
}

vector_fp SolutionArray::getComponent(const std::string& name) const
{
    if (!hasComponent(name)) {
        throw CanteraError("SolutionArray::getComponent", "no component named " + name);
    }

    vector_fp out(m_size);
    if (m_other.count(name)) {
        // auxiliary data
        auto other = m_other.at(name);
        std::copy(other->begin(), other->end(), out.begin());
        return out;
    }

    size_t ix = m_sol->thermo()->speciesIndex(name);
    if (ix == npos) {
        ix = m_sol->thermo()->nativeState()[name];
    } else {
        ix += m_stride - m_sol->thermo()->nSpecies();
    }
    for (size_t k = 0; k < m_size; ++k) {
        out[k] = m_data[k * m_stride + ix];
    }
    return out;
}

void SolutionArray::setComponent(
    const std::string& name, const vector_fp& data, bool force)
{
    if (!hasComponent(name)) {
        if (force) {
            m_other.emplace(name, std::make_shared<vector_fp>(m_size));
            auto& extra = m_other[name];
            std::copy(data.begin(), data.end(), extra->begin());
            return;
        }
        throw CanteraError("SolutionArray::setComponent", "no component named " + name);
    }
    if (data.size() != m_size) {
        throw CanteraError("SolutionArray::setComponent", "incompatible sizes");
    }

    if (m_other.count(name)) {
        // auxiliary data
        auto other = m_other[name];
        std::copy(data.begin(), data.end(), other->begin());
    }

    size_t ix = m_sol->thermo()->speciesIndex(name);
    if (ix == npos) {
        ix = m_sol->thermo()->nativeState()[name];
    } else {
        ix += m_stride - m_sol->thermo()->nSpecies();
    }
    for (size_t k = 0; k < m_size; ++k) {
        m_data[k * m_stride + ix] = data[k];
    }
}

void SolutionArray::setIndex(size_t index)
{
    if (m_size == 0) {
        throw CanteraError("SolutionArray::setIndex",
            "Unable to set index in empty SolutionArray.");
    } else if (index == npos) {
        if (m_index == npos) {
            throw CanteraError("SolutionArray::setIndex",
                "Both current and buffered indices are invalid.");
        }
        return;
    } else if (index == m_index) {
        return;
    } else if (index >= m_size) {
        throw IndexError("SolutionArray::setIndex", "entries", index, m_size - 1);
    }
    m_index = index;
    size_t nState = m_sol->thermo()->stateSize();
    m_sol->thermo()->restoreState(nState, &m_data[m_index * m_stride]);
}

vector_fp SolutionArray::getState(size_t index)
{
    setIndex(index);
    size_t nState = m_sol->thermo()->stateSize();
    vector_fp out(nState);
    m_sol->thermo()->saveState(out); // thermo contains current state
    return out;
}

std::map<std::string, double> SolutionArray::getAuxiliary(size_t index)
{
    setIndex(index);
    std::map<std::string, double> out;
    for (auto& item : m_other) {
        auto& extra = *item.second;
        out[item.first] = extra[m_index];
    }
    return out;
}

void SolutionArray::save(const std::string& fname, const std::string& id)
{
    throw CanteraError("SolutionArray::save", "Not implemented.");
}

AnyMap SolutionArray::readHeader(const std::string& fname, const std::string& id)
{
    size_t dot = fname.find_last_of(".");
    std::string extension = (dot != npos) ? toLowerCopy(fname.substr(dot + 1)) : "";
    if (extension == "h5" || extension == "hdf") {
#if CT_USE_HIGHFIVE_HDF
        return readHeader(h5::File(fname, h5::File::ReadOnly), id);
#else
        throw CanteraError("SolutionArray::readHeader",
                           "Restoring from HDF requires HighFive installation.");
#endif
    }
    if (extension == "yaml" || extension == "yml") {
        return readHeader(AnyMap::fromYamlFile(fname), id);
    }
    throw CanteraError("SolutionArray::readHeader",
                       "Unknown file extension '{}'", extension);
}

#if CT_USE_HIGHFIVE_HDF
h5::Group locateH5Group(const h5::File& file, const std::string& id)
{
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    std::string grp = tokens[0];
    if (!file.exist(grp) || file.getObjectType(grp) != h5::ObjectType::Group) {
        throw CanteraError("locateH5Group",
            "No group or solution with id '{}'", grp);
    }

    std::string path = grp;
    h5::Group sub = file.getGroup(grp);
    tokens.erase(tokens.begin());
    for (auto& grp : tokens) {
        path += "/" + grp;
        if (!sub.exist(grp) || sub.getObjectType(grp) != h5::ObjectType::Group) {
            throw CanteraError("locateH5Group",
                "No group or solution with id '{}'", path);
        }
        sub = sub.getGroup(grp);
    }
    return sub;
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
            H5Boolean value;
            attr.read(value);
            out[name] = bool(value);
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

AnyMap SolutionArray::readHeader(const h5::File& file, const std::string& id)
{
    return readH5Attributes(locateH5Group(file, id), false);
}
#endif

AnyMap SolutionArray::readHeader(const AnyMap& root, const std::string& id)
{
    throw CanteraError("SolutionArray::readHeader", "Not implemented.");
}

void SolutionArray::restore(const std::string& fname, const std::string& id)
{
    size_t dot = fname.find_last_of(".");
    std::string extension = (dot != npos) ? toLowerCopy(fname.substr(dot + 1)) : "";
    if (extension == "h5" || extension == "hdf") {
#if CT_USE_HIGHFIVE_HDF
        restore(h5::File(fname, h5::File::ReadOnly), id);
#else
        throw CanteraError("SolutionArray::restore",
                           "Restoring from HDF requires HighFive installation.");
#endif
    } else if (extension == "yaml" || extension == "yml") {
        restore(AnyMap::fromYamlFile(fname), id);
    } else {
        throw CanteraError("SolutionArray::restore",
                           "Unknown file extension '{}'", extension);
    }
}

#if CT_USE_HIGHFIVE_HDF
vector_fp readH5FloatVector(h5::DataSet data, std::string id, size_t size)
{
    if (data.getDataType().getClass() != h5::DataTypeClass::Float) {
        throw CanteraError("readH5FloatVector",
            "Type of DataSet '{}' is inconsistent; expected HDF float.", id);
    }
    if (data.getElementCount() != size) {
        throw CanteraError("readH5FloatVector",
            "Size of DataSet '{}' is inconsistent; expected {} elements but "
            "received {} elements.", id, size, data.getElementCount());
    }
    vector_fp out;
    data.read(out);
    return out;
}

std::vector<vector_fp> readH5FloatMatrix(h5::DataSet data, std::string id,
                                         size_t rows, size_t cols)
{
    if (data.getDataType().getClass() != h5::DataTypeClass::Float) {
        throw CanteraError("readH5FloatMatrix",
            "Type of DataSet '{}' is inconsistent; expected HDF float.", id);
    }
    h5::DataSpace space = data.getSpace();
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
    data.read(out);
    return out;
}

void SolutionArray::restore(const h5::File& file, const std::string& id)
{
    auto sub = locateH5Group(file, id);

    std::set<std::string> names;
    size_t nDims = npos;
    for (auto& name : sub.listObjectNames()) {
        if (sub.getObjectType(name) == h5::ObjectType::Dataset) {
            h5::DataSpace space = sub.getDataSet(name).getSpace();
            names.insert(name);
            if (space.getNumberDimensions() < nDims) {
                nDims = space.getNumberDimensions();
                m_size = space.getElementCount();
            }
        }
    }
    if (nDims != 1) {
        throw NotImplementedError("SolutionArray::restore",
            "Unable to restore SolutionArray with {} dimensions.", nDims);
    }

    initialize({});

    m_meta = readH5Attributes(sub, true);

    if (m_size == 0) {
        return;
    }

    // identify storage mode of state data
    std::string mode = "";
    const auto& nativeState = m_sol->thermo()->nativeState();
    bool usesNativeState;
    std::set<std::string> state;
    for (const auto& item : m_sol->thermo()->fullStates()) {
        bool found = true;
        usesNativeState = true;
        state.clear();
        for (size_t i = 0; i < item.size(); i++) {
            std::string name(1, item[i]);
            if (names.count(name)) {
                state.insert(name);
                usesNativeState &= nativeState.count(name);
            } else {
                found = false;
                break;
            }
        }
        if (found) {
            mode = item;
            break;
        }
    }
    if (mode == "") {
        throw CanteraError("SolutionArray::restore",
            "Data are not consistent with full state modes.");
    }

    // restore state data
    size_t nSpecies = m_sol->thermo()->nSpecies();
    size_t nState = m_sol->thermo()->stateSize();
    if (usesNativeState) {
        // native state can be written directly into data storage
        for (const auto& name : state) {
            h5::DataSet data = sub.getDataSet(name);
            if (name == "X" || name == "Y") {
                size_t offset = nativeState.find(name)->second;
                auto prop = readH5FloatMatrix(data, name, m_size, nSpecies);
                for (size_t i = 0; i < m_size; i++) {
                    std::copy(prop[i].begin(), prop[i].end(),
                              &m_data[offset + i * m_stride]);
                }
            } else {
                setComponent(name, readH5FloatVector(data, name, m_size));
            }
        }
    } else if (mode == "TPX") {
        // data format used by Python h5py export (Cantera 2.5)
        vector_fp T = readH5FloatVector(sub.getDataSet("T"), "T", m_size);
        vector_fp P = readH5FloatVector(sub.getDataSet("P"), "P", m_size);
        auto X = readH5FloatMatrix(sub.getDataSet("X"), "X", m_size, nSpecies);
        for (size_t i = 0; i < m_size; i++) {
            m_sol->thermo()->setState_TPX(T[i], P[i], X[i].data());
            m_sol->thermo()->saveState(nState, &m_data[i * m_stride]);
        }
    } else {
        throw NotImplementedError("SolutionArray::restore",
            "Import of '{}' data is not supported.", mode);
    }

    // restore other data
    for (const auto& name : names) {
        if (!state.count(name)) {
            vector_fp data = readH5FloatVector(sub.getDataSet(name), name, m_size);
            m_other.emplace(name, std::make_shared<vector_fp>(m_size));
            auto& extra = m_other[name];
            std::copy(data.begin(), data.end(), extra->begin());
        }
    }
}
#endif

void SolutionArray::restore(const AnyMap& root, const std::string& id)
{
    // locate SolutionArray based on 'id'
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    const AnyMap* ptr = &root; // use raw pointer to avoid copying
    std::string path = "";
    for (auto& field : tokens) {
        path += "/" + field;
        const AnyMap& sub = *ptr;
        if (!sub.hasKey(field) || !sub[field].is<AnyMap>()) {
            throw CanteraError("SolutionArray::restore",
                "No field or solution with id '{}'", path);
        }
        ptr = &sub[field].as<AnyMap>(); // AnyMap lacks 'operator=' for const AnyMap
    }
    const AnyMap& sub = *ptr;

    // set size and initialize
    m_size = sub.getInt("points", 0);
    if (!sub.hasKey("T") && !sub.hasKey("temperature")) {
        // overwrite size - Sim1D erroneously assigns '1' (Cantera 2.6)
        m_size = 0;
    }
    initialize({});

    // restore data
    std::set<std::string> exclude = {"points", "X", "Y"};
    if (m_size == 0) {
        // no data points
    } else if (m_size == 1) {
        // single data point
        double T = sub["temperature"].asDouble();
        double P = sub.getDouble("pressure", OneAtm); // missing - Sim1D (Cantera 2.6)
        std::set<std::string> props = {"temperature", "pressure"};
        exclude.insert(props.begin(), props.end());
        if (sub.hasKey("mass-fractions")) {
            auto Y = sub["mass-fractions"].asMap<double>();
            m_sol->thermo()->setState_TPY(T, P, Y);
            exclude.insert("mass-fractions");
        } else if (sub.hasKey("coverages")) {
            m_sol->thermo()->setState_TP(T, P);
            auto cov = sub["coverages"].asMap<double>();
            exclude.insert("coverages");
            auto surf = std::dynamic_pointer_cast<SurfPhase>(m_sol->thermo());
            if (!surf) {
                throw CanteraError("SolutionArray::restore",
                    "Restoring of coverages requires surface phase");
            }
            surf->setCoveragesByName(cov);
        } else {
            throw NotImplementedError("SolutionArray::restore",
                "Unknown YAML serialization format.");
        }
        for (const auto& prop : m_sol->thermo()->nativeState()) {
            if (prop.first == "T") {
                m_data[prop.second] = m_sol->thermo()->temperature();
            } else if (prop.first == "D") {
                m_data[prop.second] = m_sol->thermo()->density();
            } else if (prop.first == "P") {
                m_data[prop.second] = m_sol->thermo()->pressure();
            } else if (prop.first == "Y") {
                m_sol->thermo()->getMassFractions(&m_data[prop.second]);
            } else if (prop.first == "X") {
                m_sol->thermo()->getMoleFractions(&m_data[prop.second]);
            } else {
                throw NotImplementedError("SolutionArray::restore",
                    "Unable to restore property '{}'.", prop.first);
            }
        }
    } else {
        // multiple data points
        const auto& nativeState = m_sol->thermo()->nativeState();
        for (const auto& item : sub) {
            const std::string& name = item.first;
            const AnyValue& value = item.second;
            if (value.is<std::vector<double>>()) {
                const vector_fp& data = value.as<std::vector<double>>();
                if (data.size() == m_size) {
                    setComponent(name, data, true);
                    exclude.insert(item.first);
                }
            }
        }

        // check that state data are complete
        std::set<std::string> props = {};
        std::set<std::string> missingProps = {};
        for (const auto& item : nativeState) {
            if (exclude.count(item.first)) {
                props.insert(item.first);
            } else {
                missingProps.insert(item.first);
            }
        }

        std::set<std::string> TY = {"T", "Y"};
        if (props == TY && missingProps.count("D") && sub.hasKey("pressure")) {
            // missing "D" - Sim1D (Cantera 2.6)
            double P = sub["pressure"].asDouble();
            const size_t offset_T = nativeState.find("T")->second;
            const size_t offset_D = nativeState.find("D")->second;
            const size_t offset_Y = nativeState.find("Y")->second;
            for (size_t i = 0; i < m_size; i++) {
                double T = m_data[offset_T + i * m_stride];
                m_sol->thermo()->setState_TPY(T, P, &m_data[offset_Y + i * m_stride]);
                m_data[offset_D + i * m_stride] = m_sol->thermo()->density();
            }
        } else if (missingProps.size()) {
            throw CanteraError("SolutionArray::restore",
                "Incomplete state information.");
        }
    }

    // add meta data
    for (const auto& item : sub) {
        if (!exclude.count(item.first)) {
            m_meta[item.first] = item.second;
        }
    }
}

}
