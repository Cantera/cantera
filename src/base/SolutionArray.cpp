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
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;
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
    size_t count = 0;
    for (auto& key : extra) {
        m_extra.emplace(key, count);
        count++;
    }

    m_offsets = m_sol->thermo()->nativeState();
    m_stride = m_sol->thermo()->stateSize();
    m_work.reset(new vector_fp(m_size * m_stride, 0.));
    m_data = m_work->data();
    m_managed = false;
    for (auto& key : extra) {
        m_other.emplace(key, std::make_shared<vector_fp>(m_size));
    }
}

void SolutionArray::initialize(
    double* data,
    size_t size,
    size_t stride,
    const std::vector<std::pair<std::string, size_t>>& offsets)
{
    // check that offsets match order of native thermodynamic state properties
    std::map<size_t, std::string> flipped;
    for (const auto& item : m_sol->thermo()->nativeState()) {
        // flipped map will be sorted by native property offset within state
        flipped.emplace(item.second, item.first);
    }
    std::map<std::string, int> mapped; // searchable offset map
    for (const auto& item : offsets) {
        mapped.emplace(item.first, (int)(item.second));
    }
    std::string key0 = flipped.at(0);
    for (auto& prop : flipped) {
        if (!mapped.count(prop.second)) {
            throw CanteraError("SolutionArray::initialize",
                "Native property '{}' not found in offset mapping.", prop.second);
        }
        int diffOffset = mapped.at(prop.second) - mapped.at(key0);
        if (diffOffset != (int)(prop.first)) {
            throw CanteraError("SolutionArray::initialize",
                "Offset for property '{}' is incompatible with order of native state "
                "properties", prop.second);
        }
    }

    // assign managed memory
    m_work.reset();
    m_data = data;
    m_managed = true;
    m_size = size;
    m_stride = stride;

    size_t count = 0;
    for (auto& item : offsets) {
        auto& key = item.first;
        if (item.second != npos) {
            m_offsets[key] = item.second;
        } else {
            m_other.emplace(key, std::make_shared<vector_fp>(m_size));
            m_extra.emplace(key, count);
            count++;
        }
    }
}

void SolutionArray::save(const std::string& fname, const std::string& id)
{
    throw CanteraError("SolutionArray::save", "Not implemented.");
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
void SolutionArray::restore(const h5::File& file, const std::string& id)
{
    std::vector<std::string> tokens;
    tokenizePath(id, tokens);
    std::string grp = tokens[0];
    if (!file.exist(grp) || file.getObjectType(grp) != h5::ObjectType::Group) {
        throw CanteraError("SolutionArray::restore",
                           "No group or solution with id '{}'", grp);
    }

    std::string path = grp;
    h5::Group sub = file.getGroup(grp);
    tokens.erase(tokens.begin());
    for (auto& grp : tokens) {
        path += "/" + grp;
        if (!sub.exist(grp) || sub.getObjectType(grp) != h5::ObjectType::Group) {
            throw CanteraError("SolutionArray::restore",
                               "No group or solution with id '{}'", path);
        }
        sub = sub.getGroup(grp);
    }

    std::vector<std::string> names;
    size_t nDims = npos;
    for (auto& name : sub.listObjectNames()) {
        if (sub.getObjectType(name) == h5::ObjectType::Dataset) {
            h5::DataSpace space = sub.getDataSet(name).getSpace();
            names.push_back(name);
            if (space.getNumberDimensions() < nDims) {
                nDims = space.getNumberDimensions();
                m_size = space.getElementCount();
            }
        }
    }

    // @todo: restore data
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
            size_t offset = npos;
            if (value.is<std::vector<double>>()) {
                const vector_fp& data = value.as<std::vector<double>>();
                size_t species = m_sol->thermo()->speciesIndex(name);
                if (data.size() != m_size) {
                    // meta data
                    continue;
                } else if (species != npos) {
                    // species
                    if (nativeState.count("X")) {
                        offset = nativeState.find("X")->second + species;
                    } else if (nativeState.count("Y")) {
                        offset = nativeState.find("Y")->second + species;
                    }
                } else if (nativeState.count(name)) {
                    // property
                    offset = nativeState.find(name)->second;
                } else {
                    // extra
                    m_other.emplace(name, std::make_shared<vector_fp>(m_size));
                    auto& extra = m_other[name];
                    std::copy(data.begin(), data.end(), extra->begin());
                }

                if (offset != npos) {
                    for (size_t i = 0; i < m_size; i++) {
                        m_data[offset + i * m_stride] = data[i];
                    }
                }
                exclude.insert(item.first);
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
