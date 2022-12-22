/**
 * @file SolutionArray.cpp
 *      Definition file for class SolutionArray.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/SolutionArray.h"
#include "cantera/base/Solution.h"
#include "cantera/base/Storage.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include <boost/algorithm/string/predicate.hpp>
#include <set>
#include <fstream>


namespace ba = boost::algorithm;

namespace Cantera
{

const std::map<std::string, std::string> aliasMap = {
    {"T", "temperature"},
    {"P", "pressure"},
    {"D", "density"},
    {"Y", "mass-fractions"},
    {"X", "mole-fractions"},
    {"C", "coverages"},
    {"U", "specific-internal-energy"},
    {"V", "specific-volume"},
    {"H", "specific-enthalpy"},
    {"S", "specific-entropy"},
    {"Q", "vapor-fraction"},
};

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
    m_stride = m_sol->thermo()->stateSize();
    m_data.resize(m_size * m_stride, 0.);
}

shared_ptr<ThermoPhase> SolutionArray::thermo()
{
    return m_sol->thermo();
}

bool SolutionArray::hasComponent(const std::string& name) const
{
    if (m_extra.count(name)) {
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
    if (m_extra.count(name)) {
        // auxiliary data
        out = m_extra.at(name);
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
            m_extra.emplace(name, data);
            return;
        }
        throw CanteraError("SolutionArray::setComponent", "no component named " + name);
    }
    if (data.size() != m_size) {
        throw CanteraError("SolutionArray::setComponent", "incompatible sizes");
    }

    if (m_extra.count(name)) {
        // auxiliary data
        m_extra[name] = data;
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

void SolutionArray::setIndex(size_t index, bool restore)
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
    if (restore) {
        size_t nState = m_sol->thermo()->stateSize();
        m_sol->thermo()->restoreState(nState, &m_data[m_index * m_stride]);
    }
}

vector_fp SolutionArray::getState(size_t index)
{
    setIndex(index);
    size_t nState = m_sol->thermo()->stateSize();
    vector_fp out(nState);
    m_sol->thermo()->saveState(out); // thermo contains current state
    return out;
}

void SolutionArray::setState(const vector_fp& data, size_t index)
{
    setIndex(index, false);
    m_sol->thermo()->restoreState(data);
    size_t nState = m_sol->thermo()->stateSize();
    m_sol->thermo()->saveState(nState, &m_data[m_index * m_stride]);
}

std::map<std::string, double> SolutionArray::getAuxiliary(size_t index)
{
    setIndex(index);
    std::map<std::string, double> out;
    for (auto& item : m_extra) {
        out[item.first] = item.second[m_index];
    }
    return out;
}

AnyMap preamble(const std::string& desc)
{
    AnyMap data;
    data["description"] = desc;
    data["generator"] = "Cantera SolutionArray";
    data["cantera-version"] = CANTERA_VERSION;
    data["git-commit"] = gitCommit();

    // Add a timestamp indicating the current time
    time_t aclock;
    ::time(&aclock); // Get time in seconds
    struct tm* newtime = localtime(&aclock); // Convert time to struct tm form
    data["date"] = stripnonprint(asctime(newtime));

    // Force metadata fields to the top of the file
    data["description"].setLoc(-6, 0);
    data["generator"].setLoc(-5, 0);
    data["cantera-version"].setLoc(-4, 0);
    data["git-commit"].setLoc(-3, 0);
    data["date"].setLoc(-2, 0);

    return data;
}

void SolutionArray::writeHeader(const std::string& fname, const std::string& id,
                                const std::string& desc)
{
    Storage file(fname, true);
    file.checkGroup(id);
    file.writeAttributes(id, preamble(desc));
}

void SolutionArray::writeHeader(AnyMap& root, const std::string& id,
                                const std::string& desc)
{
    root[id] = preamble(desc);
}

void SolutionArray::writeEntry(const std::string& fname, const std::string& id,
                               int compression)
{
    Storage file(fname, true);
    if (compression) {
        file.setCompressionLevel(compression);
    }
    file.checkGroup(id);
    file.writeAttributes(id, m_meta);
    if (!m_size) {
        return;
    }

    const auto& nativeState = m_sol->thermo()->nativeState();
    size_t nSpecies = m_sol->thermo()->nSpecies();
    for (auto& state : nativeState) {
        std::string name = state.first;
        if (name == "X" || name == "Y") {
            size_t offset = state.second;
            std::vector<vector_fp> prop;
            for (size_t i = 0; i < m_size; i++) {
                size_t first = offset + i * m_stride;
                prop.push_back(vector_fp(&m_data[first], &m_data[first + nSpecies]));
            }
            file.writeMatrix(id, name, prop);
        } else {
            auto data = getComponent(name);
            file.writeVector(id, name, data);
        }
    }

    for (auto& extra : m_extra) {
        file.writeVector(id, extra.first, extra.second);
    }
}

AnyMap& openField(AnyMap& root, const std::string& id)
{
    // locate field based on 'id'
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    AnyMap* ptr = &root; // use raw pointer to avoid copying
    std::string path = "";
    for (auto& field : tokens) {
        path += "/" + field;
        AnyMap& sub = *ptr;
        if (sub.hasKey(field) && !sub[field].is<AnyMap>()) {
            throw CanteraError("openField",
                "Encountered invalid existing field '{}'", path);
        } else if (!sub.hasKey(field)) {
            sub[field] = AnyMap();
        }
        ptr = &sub[field].as<AnyMap>();
    }
    return *ptr;
}

void SolutionArray::writeEntry(AnyMap& root, const std::string& id)
{
    AnyMap& data = openField(root, id);
    bool preexisting = !data.empty();
    data["points"] = int(m_size);
    data.update(m_meta);

    for (auto& extra : m_extra) {
        data[extra.first] = extra.second;
    }

    auto phase = m_sol->thermo();
    if (m_size == 1) {
        setIndex(0);
        data["temperature"] = phase->temperature();
        data["pressure"] = phase->pressure();
        auto surf = std::dynamic_pointer_cast<SurfPhase>(phase);
        auto nSpecies = phase->nSpecies();
        vector_fp values(nSpecies);
        if (surf) {
            surf->getCoverages(&values[0]);
        } else {
            phase->getMassFractions(&values[0]);
        }
        AnyMap items;
        for (size_t k = 0; k < nSpecies; k++) {
            if (values[k] != 0.0) {
                items[phase->speciesName(k)] = values[k];
            }
        }
        if (surf) {
            data["coverages"] = std::move(items);
        } else {
            data["mass-fractions"] = std::move(items);
        }
    } else if (m_size > 1) {
        const auto& nativeState = phase->nativeState();
        for (auto& state : nativeState) {
            std::string name = state.first;
            if (name == "X" || name == "Y") {
                data["basis"] = name == "X" ? "mole" : "mass";
                for (auto& name : phase->speciesNames()) {
                    data[name] = getComponent(name);
                }
            } else {
                data[name] = getComponent(name);
            }
        }
    }

    // If this is not replacing an existing solution, put it at the end
    if (!preexisting) {
        data.setLoc(INT_MAX, 0);
    }
}

void SolutionArray::save(const std::string& fname, const std::string& id,
                         const std::string& desc, int compression)
{
    size_t dot = fname.find_last_of(".");
    std::string extension = (dot != npos) ? toLowerCopy(fname.substr(dot + 1)) : "";
    if (extension == "h5" || extension == "hdf"  || extension == "hdf5") {
        writeHeader(fname, id, desc);
        writeEntry(fname, id, compression);
        return;
    }
    if (extension == "yaml" || extension == "yml") {
        // Check for an existing file and load it if present
        AnyMap data;
        if (std::ifstream(fname).good()) {
            data = AnyMap::fromYamlFile(fname);
        }
        writeEntry(data, id);
        writeHeader(data, id, desc);

        // Write the output file and remove the now-outdated cached file
        std::ofstream out(fname);
        out << data.toYamlString();
        AnyMap::clearCachedFile(fname);
        return;
    }
    throw CanteraError("SolutionArray::writeHeader",
                       "Unknown file extension '{}'", extension);
}

AnyMap SolutionArray::readHeader(const std::string& fname, const std::string& id)
{
    Storage file(fname, false);
    file.checkGroup(id);
    return file.readAttributes(id, false);
}

const AnyMap& locateField(const AnyMap& root, const std::string& id)
{
    // locate field based on 'id'
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
    return *ptr;
}

AnyMap SolutionArray::readHeader(const AnyMap& root, const std::string& id)
{
    auto sub = locateField(root, id);
    AnyMap header;
    for (const auto& item : sub) {
        if (!sub[item.first].is<AnyMap>()) {
            header[item.first] = item.second;
        }
    }
    return header;
}

AnyMap SolutionArray::restore(const std::string& fname, const std::string& id)
{
    size_t dot = fname.find_last_of(".");
    std::string extension = (dot != npos) ? toLowerCopy(fname.substr(dot + 1)) : "";
    if (extension == "h5" || extension == "hdf"  || extension == "hdf5") {
        readEntry(fname, id);
        return readHeader(fname, id);
    }
    if (extension == "yaml" || extension == "yml") {
        const AnyMap& root = AnyMap::fromYamlFile(fname);
        readEntry(root, id);
        return readHeader(root, id);
    }
    throw CanteraError("SolutionArray::restore",
                        "Unknown file extension '{}'", extension);
}

std::string SolutionArray::detectMode(const std::set<std::string>& names, bool native)
{
    // check set of available names against state acronyms defined by Phase::fullStates
    std::string mode = "";
    const auto& nativeState = m_sol->thermo()->nativeState();
    bool usesNativeState;
    auto surf = std::dynamic_pointer_cast<SurfPhase>(m_sol->thermo());
    for (const auto& item : m_sol->thermo()->fullStates()) {
        bool found = true;
        std::string name;
        usesNativeState = true;
        for (size_t i = 0; i < item.size(); i++) {
            // pick i-th letter from "full" state acronym
            name = std::string(1, item[i]);
            if (surf && (name == "X" || name == "Y")) {
                // override native state to enable detection of surface phases
                name = "C";
                usesNativeState = false;
                break;
            }
            if (names.count(name)) {
                // property is stored using letter acronym
                usesNativeState &= nativeState.count(name) > 0;
            } else if (aliasMap.count(name) && names.count(aliasMap.at(name))) {
                // property is stored using property name
                usesNativeState &= nativeState.count(name) > 0;
            } else {
                found = false;
                break;
            }
        }
        if (found) {
            mode = (name == "C") ? item.substr(0, 2) + "C" : item;
            break;
        }
    }
    if (usesNativeState && native) {
        return "native";
    }
    return mode;
}

std::set<std::string> SolutionArray::stateProperties(
    const std::string& mode, bool alias)
{
    std::set<std::string> states;
    if (mode == "native") {
        for (const auto& item : m_sol->thermo()->nativeState()) {
            states.insert(alias ? aliasMap.at(item.first) : item.first);
        }
    } else {
        for (const auto& m : mode) {
            const std::string name = std::string(1, m);
            states.insert(alias ? aliasMap.at(name) : name);
        }
    }

    return states;
}

void SolutionArray::readEntry(const std::string& fname, const std::string& id)
{
    Storage file(fname, false);
    file.checkGroup(id);
    m_meta = file.readAttributes(id, true);

    auto contents = file.contents(id);
    m_size = contents.first;
    m_data.resize(m_size * m_stride, 0.);
    std::set<std::string> names = contents.second;

    if (m_size == 0) {
        return;
    }

    // determine storage mode of state data
    std::string mode = detectMode(names);
    std::set<std::string> states = stateProperties(mode);
    if (states.count("C")) {
        if (names.count("X")) {
            states.erase("C");
            states.insert("X");
            mode = "TPX";
        } else if (names.count("Y")) {
            states.erase("C");
            states.insert("Y");
            mode = "TPY";
        }
    }

    // restore state data
    size_t nSpecies = m_sol->thermo()->nSpecies();
    size_t nState = m_sol->thermo()->stateSize();
    const auto& nativeStates = m_sol->thermo()->nativeState();
    if (mode == "native") {
        // native state can be written directly into data storage
        for (const auto& item : nativeStates) {
            std::string name = item.first;
            if (name == "X" || name == "Y") {
                size_t offset = item.second;
                auto prop = file.readMatrix(id, name, m_size, nSpecies);
                for (size_t i = 0; i < m_size; i++) {
                    std::copy(prop[i].begin(), prop[i].end(),
                              &m_data[offset + i * m_stride]);
                }
            } else {
                setComponent(name, file.readVector(id, name, m_size));
            }
        }
    } else if (mode == "TPX") {
        // data format used by Python h5py export (Cantera 2.5)
        vector_fp T = file.readVector(id, "T", m_size);
        vector_fp P = file.readVector(id, "P", m_size);
        auto X = file.readMatrix(id, "X", m_size, nSpecies);
        for (size_t i = 0; i < m_size; i++) {
            m_sol->thermo()->setState_TPX(T[i], P[i], X[i].data());
            m_sol->thermo()->saveState(nState, &m_data[i * m_stride]);
        }
    } else if (mode == "TPY") {
        vector_fp T = file.readVector(id, "T", m_size);
        vector_fp P = file.readVector(id, "P", m_size);
        auto Y = file.readMatrix(id, "Y", m_size, nSpecies);
        for (size_t i = 0; i < m_size; i++) {
            m_sol->thermo()->setState_TPY(T[i], P[i], Y[i].data());
            m_sol->thermo()->saveState(nState, &m_data[i * m_stride]);
        }
    } else if (mode == "") {
        throw CanteraError("SolutionArray::restore",
            "Data are not consistent with full state modes.");
    } else {
        throw NotImplementedError("SolutionArray::restore",
            "Import of '{}' data is not supported.", mode);
    }

    // restore remaining data
    for (const auto& name : names) {
        if (!states.count(name)) {
            vector_fp data = file.readVector(id, name, m_size);
            m_extra.emplace(name, data);
        }
    }
}

void SolutionArray::readEntry(const AnyMap& root, const std::string& id)
{
    auto sub = locateField(root, id);

    // set size and initialize
    m_size = sub.getInt("points", 0);
    if (!sub.hasKey("T") && !sub.hasKey("temperature")) {
        // overwrite size - Sim1D erroneously assigns '1' (Cantera 2.6)
        m_size = 0;
    }
    m_data.resize(m_size * m_stride, 0.);

    // restore data
    std::set<std::string> exclude = {"points", "X", "Y"};
    std::set<std::string> names = sub.keys();
    size_t nState = m_sol->thermo()->stateSize();
    if (m_size == 0) {
        // no data points
    } else if (m_size == 1) {
        // single data point
        std::string mode = detectMode(names, false);
        if (mode == "") {
            // missing property - Sim1D (Cantera 2.6)
            names.insert("pressure");
            mode = detectMode(names, false);
        }
        if (mode == "TPY") {
            // single data point uses long names
            double T = sub["temperature"].asDouble();
            double P = sub.getDouble("pressure", OneAtm); // missing - Sim1D (Cantera 2.6)
            auto Y = sub["mass-fractions"].asMap<double>();
            m_sol->thermo()->setState_TPY(T, P, Y);
        } else if (mode == "TPC") {
            auto surf = std::dynamic_pointer_cast<SurfPhase>(m_sol->thermo());
            if (!surf) {
                throw CanteraError("SolutionArray::restore",
                    "Restoring of coverages requires surface phase");
            }
            double T = sub["temperature"].asDouble();
            double P = sub.getDouble("pressure", OneAtm); // missing - Sim1D (Cantera 2.6)
            m_sol->thermo()->setState_TP(T, P);
            auto cov = sub["coverages"].asMap<double>();
            surf->setCoveragesByName(cov);
        } else if (mode == "") {
            throw CanteraError("SolutionArray::restore",
                "Data are not consistent with full state modes.");
        } else {
            throw NotImplementedError("SolutionArray::restore",
                "Import of '{}' data is not supported.", mode);
        }
        m_sol->thermo()->saveState(nState, m_data.data());
        auto props = stateProperties(mode, true);
        exclude.insert(props.begin(), props.end());
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
                "Incomplete state information: missing '{}'",
                ba::join(missingProps, "', '"));
        }
    }

    // add meta data
    for (const auto& item : sub) {
        if (!exclude.count(item.first)) {
            m_meta[item.first] = item.second;
        }
    }
    m_meta.erase("points");
}

}
