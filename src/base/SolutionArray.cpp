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
#include "cantera/base/utilities.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>


namespace ba = boost::algorithm;

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

namespace Cantera
{

SolutionArray::SolutionArray(const shared_ptr<Solution>& sol,
                             int size, const AnyMap& meta)
    : m_sol(sol)
    , m_size(size)
    , m_dataSize(size)
    , m_meta(meta)
{
    if (!m_sol || !m_sol->thermo()) {
        throw CanteraError("SolutionArray::SolutionArray",
            "Unable to create SolutionArray from invalid Solution object.");
    }
    m_stride = m_sol->thermo()->stateSize();
    m_sol->thermo()->addSpeciesLock();
    m_data = make_shared<vector<double>>(m_dataSize * m_stride, 0.);
    m_extra = make_shared<map<string, AnyValue>>();
    m_order = make_shared<map<int, string>>();
    for (size_t i = 0; i < m_dataSize; ++i) {
        m_active.push_back(static_cast<int>(i));
    }
    reset();
    m_apiShape.resize(1);
    m_apiShape[0] = static_cast<long>(m_dataSize);
}

SolutionArray::SolutionArray(const SolutionArray& other,
                             const vector<int>& selected)
    : m_sol(other.m_sol)
    , m_size(selected.size())
    , m_dataSize(other.m_data->size())
    , m_stride(other.m_stride)
    , m_data(other.m_data)
    , m_extra(other.m_extra)
    , m_order(other.m_order)
    , m_shared(true)
{
    m_sol->thermo()->addSpeciesLock();
    if (!other.m_shared) {
        // direct slicing is possible
        m_active = selected;
    } else {
        // slicing a previously sliced SolutionArray
        m_active.clear();
        m_active.reserve(selected.size());
        for (auto loc : selected) {
            m_active.push_back(other.m_active.at(loc));
        }
    }
    for (auto loc : m_active) {
        if (loc < 0 || loc >= (int)m_dataSize) {
            IndexError("SolutionArray::SolutionArray", "indices", loc, m_dataSize);
        }
    }
    set<int> unique(selected.begin(), selected.end());
    if (unique.size() < selected.size()) {
        throw CanteraError("SolutionArray::SolutionArray", "Indices must be unique.");
    }
}

SolutionArray::~SolutionArray()
{
    m_sol->thermo()->removeSpeciesLock();
}

namespace { // restrict scope of helper functions to local translation unit

template<class T>
void resetSingle(AnyValue& extra, const vector<int>& slice);

template<class T>
AnyValue getSingle(const AnyValue& extra, const vector<int>& slice);

template<class T>
void setSingle(AnyValue& extra, const AnyValue& data, const vector<int>& slice);

template<class T>
void resizeSingle(AnyValue& extra, size_t size, const AnyValue& value);

template<class T>
void resetMulti(AnyValue& extra, const vector<int>& slice);

template<class T>
AnyValue getMulti(const AnyValue& extra, const vector<int>& slice);

template<class T>
void setMulti(AnyValue& extra, const AnyValue& data, const vector<int>& slice);

template<class T>
void resizeMulti(AnyValue& extra, size_t size, const AnyValue& value);

template<class T>
void setAuxiliarySingle(size_t loc, AnyValue& extra, const AnyValue& value);

template<class T>
void setAuxiliaryMulti(size_t loc, AnyValue& extra, const AnyValue& data);

template<class T>
void setScalar(AnyValue& extra, const AnyValue& data, const vector<int>& slice);

} // end unnamed namespace

void SolutionArray::reset()
{
    size_t nState = m_sol->thermo()->stateSize();
    vector<double> state(nState);
    m_sol->thermo()->saveState(state); // thermo contains current state
    for (size_t k = 0; k < m_size; ++k) {
        std::copy(state.begin(), state.end(), m_data->data() + m_active[k] * m_stride);
    }
    for (auto& [key, extra] : *m_extra) {
        if (extra.is<void>()) {
            // cannot reset placeholder (uninitialized component)
        } else if (extra.isVector<double>()) {
            resetSingle<double>(extra, m_active);
        } else if (extra.isVector<long int>()) {
            resetSingle<long int>(extra, m_active);
        } else if (extra.isVector<string>()) {
            resetSingle<string>(extra, m_active);
        } else if (extra.isVector<vector<double>>()) {
            resetMulti<double>(extra, m_active);
        } else if (extra.isVector<vector<long int>>()) {
            resetMulti<long int>(extra, m_active);
        } else if (extra.isVector<vector<string>>()) {
            resetMulti<string>(extra, m_active);
        } else {
            throw NotImplementedError("SolutionArray::reset",
                "Unable to reset component '{}' with type '{}'.",
                key, extra.type_str());
        }
    }
}

void SolutionArray::resize(int size)
{
    if (apiNdim() > 1) {
        throw CanteraError("SolutionArray::resize",
            "Resize is ambiguous for multi-dimensional arrays; use setApiShape "
            "instead.");
    }
    if (m_data.use_count() > 1) {
        throw CanteraError("SolutionArray::resize",
            "Unable to resize as data are shared by multiple objects.");
    }
    _resize(static_cast<size_t>(size));
    m_apiShape[0] = static_cast<long>(size);
}

void SolutionArray::setApiShape(const vector<long int>& shape)
{
    size_t size = 1;
    for (auto dim : shape) {
        size *= dim;
    }
    if (m_shared && size != m_size) {
        throw CanteraError("SolutionArray::setApiShape",
            "Unable to set shape of shared data as sizes are inconsistent:\n"
            "active size is {} but shape implies {}.", m_size, size);
    }
    if (!m_shared && size != m_dataSize) {
        if (m_data.use_count() > 1) {
            throw CanteraError("SolutionArray::setApiShape",
                "Unable to set shape as data are shared by multiple objects.");
        }
        _resize(size);
    }
    m_apiShape = shape;
}

void SolutionArray::_resize(size_t size)
{
    m_size = size;
    m_dataSize = size;
    m_data->resize(m_dataSize * m_stride, 0.);
    for (auto& [key, data] : *m_extra) {
        _resizeExtra(key);
    }
    m_active.clear();
    for (size_t i = 0; i < m_dataSize; ++i) {
        m_active.push_back(static_cast<int>(i));
    }
}

namespace { // restrict scope of helper functions to local translation unit

vector<string> doubleColumn(string name, const vector<double>& comp,
                            int rows, int width)
{
    // extract data for processing
    vector<double> data;
    vector<string> raw;
    string notation = fmt::format("{{:{}.{}g}}", width, (width - 1) / 2);
    int csize = static_cast<int>(comp.size());
    int dots = csize + 1;
    if (csize <= rows) {
        for (const auto& val : comp) {
            data.push_back(val);
            raw.push_back(boost::trim_copy(fmt::format(notation, val)));
        }
    } else {
        dots = (rows + 1) / 2;
        for (int row = 0; row < dots; row++) {
            data.push_back(comp[row]);
            raw.push_back(boost::trim_copy(fmt::format(notation, comp[row])));
        }
        for (int row = csize - rows / 2; row < csize; row++) {
            data.push_back(comp[row]);
            raw.push_back(boost::trim_copy(fmt::format(notation, comp[row])));
        }
    }

    // determine notation; all entries use identical formatting
    bool isFloat = false;
    bool isScientific = false;
    size_t head = 0;
    size_t tail = 0;
    size_t exp = 0;
    for (const auto& repr : raw) {
        string name = repr;
        if (name[0] == '-') {
            // leading negative sign is not considered
            name = name.substr(1);
        }
        if (name.find('e') != string::npos) {
            // scientific notation
            if (!isScientific) {
                head = 1;
                tail = name.find('e') - name.find('.');
                exp = 4; // size of exponent
            } else {
                tail = std::max(tail, name.find('e') - name.find('.'));
            }
            isFloat = true;
            isScientific = true;
        } else if (name.find('.') != string::npos) {
            // floating point notation
            isFloat = true;
            if (!isScientific) {
                head = std::max(head, name.find('.'));
                tail = std::max(tail, name.size() - name.find('.'));
            }
        } else {
            head = std::max(head, name.size());
        }
    }
    size_t maxLen = std::max((size_t)4, head + tail + exp + isFloat + 1);
    size_t over = std::max(0, (int)name.size() - (int)maxLen);
    if (isScientific) {
        // at least one entry has scientific notation
        notation = fmt::format(" {{:>{}.{}e}}", over + maxLen, tail);
    } else if (isFloat) {
        // at least one entry is a floating point
        notation = fmt::format(" {{:>{}.{}f}}", over + maxLen, tail);
    } else {
        // all entries are integers
        notation = fmt::format(" {{:>{}.0f}}", over + maxLen);
    }
    maxLen = fmt::format(notation, 0.).size();

    // assemble output
    string section = fmt::format("{{:>{}}}", maxLen);
    vector<string> col = {fmt::format(section, name)};
    int count = 0;
    for (const auto& val : data) {
        col.push_back(fmt::format(notation, val));
        count++;
        if (count == dots) {
            col.push_back(fmt::format(section, "..."));
        }
    }
    return col;
}

vector<string> integerColumn(string name, const vector<long int>& comp,
                             int rows, int width)
{
    // extract data for processing
    vector<long int> data;
    string notation = fmt::format("{{:{}}}", width);
    size_t maxLen = 2; // minimum column width is 2
    int csize = static_cast<int>(comp.size());
    int dots = csize + 1;
    if (csize <= rows) {
        for (const auto& val : comp) {
            data.push_back(val);
            string formatted = boost::trim_copy(fmt::format(notation, val));
            if (formatted[0] == '-') {
                formatted = formatted.substr(1);
            }
            maxLen = std::max(maxLen, formatted.size());
        }
    } else {
        dots = (rows + 1) / 2;
        for (int row = 0; row < dots; row++) {
            data.push_back(comp[row]);
            string formatted = boost::trim_copy(fmt::format(notation, comp[row]));
            if (formatted[0] == '-') {
                formatted = formatted.substr(1);
            }
            maxLen = std::max(maxLen, formatted.size());
        }
        for (int row = csize - rows / 2; row < csize; row++) {
            data.push_back(comp[row]);
            string formatted = boost::trim_copy(fmt::format(notation, comp[row]));
            if (formatted[0] == '-') {
                formatted = formatted.substr(1);
            }
            maxLen = std::max(maxLen, formatted.size());
        }
    }

    if (name == "") {
        // index column
        notation = fmt::format("{{:<{}}}", maxLen);
    } else {
        // regular column
        maxLen = std::max(maxLen, name.size());
        notation = fmt::format(" {{:>{}}}", maxLen + 1);
    }

    // assemble output
    vector<string> col = {fmt::format(notation, name)};
    int count = 0;
    for (const auto& val : data) {
        col.push_back(fmt::format(notation, val));
        count++;
        if (count == dots) {
            col.push_back(fmt::format(notation, ".."));
        }
    }
    return col;
}

vector<string> stringColumn(string name, const vector<string>& comp,
                            int rows, int width)
{
    // extract data for processing
    vector<string> data;
    string notation = fmt::format("{{:{}}}", width);
    size_t maxLen = 3; // minimum column width is 3
    int csize = static_cast<int>(comp.size());
    int dots = csize + 1;
    if (csize <= rows) {
        for (const auto& val : comp) {
            data.push_back(val);
            maxLen = std::max(maxLen,
                boost::trim_copy(fmt::format(notation, val)).size());
        }
    } else {
        dots = (rows + 1) / 2;
        for (int row = 0; row < dots; row++) {
            data.push_back(comp[row]);
            maxLen = std::max(maxLen,
                boost::trim_copy(fmt::format(notation, comp[row])).size());
        }
        for (int row = csize - rows / 2; row < csize; row++) {
            data.push_back(comp[row]);
            maxLen = std::max(maxLen,
                boost::trim_copy(fmt::format(notation, comp[row])).size());
        }
    }

    // assemble output
    notation = fmt::format("  {{:>{}}}", maxLen);
    vector<string> col = {fmt::format(notation, name)};
    int count = 0;
    for (const auto& val : data) {
        col.push_back(fmt::format(notation, val));
        count++;
        if (count == dots) {
            col.push_back(fmt::format(notation, "..."));
        }
    }
    return col;
}

vector<string> formatColumn(string name, const AnyValue& comp, int rows, int width)
{
    if (comp.isVector<double>()) {
        return doubleColumn(name, comp.asVector<double>(), rows, width);
    }
    if (comp.isVector<long int>()) {
        return integerColumn(name, comp.asVector<long int>(), rows, width);
    }
    if (comp.isVector<string>()) {
        return stringColumn(name, comp.asVector<string>(), rows, width);
    }

    // create alternative representation
    string repr;
    int size;
    if (comp.isVector<vector<double>>()) {
        repr = "[ <double> ]";
        size = len(comp.asVector<vector<double>>());
    } else if (comp.isVector<vector<long int>>()) {
        repr = "[ <long int> ]";
        size = len(comp.asVector<vector<long int>>());
    } else if (comp.isVector<vector<string>>()) {
        repr = "[ <string> ]";
        size = len(comp.asVector<vector<string>>());
    } else {
        throw CanteraError(
            "formatColumn", "Encountered invalid data for component '{}'.", name);
    }
    size_t maxLen = std::max(repr.size(), name.size());

    // assemble output
    string notation = fmt::format("  {{:>{}}}", maxLen);
    repr = fmt::format(notation, repr);
    vector<string> col = {fmt::format(notation, name)};
    if (size <= rows) {
        for (int row = 0; row < size; row++) {
            col.push_back(repr);
        }
    } else {
        int dots = (rows + 1) / 2;
        for (int row = 0; row < dots; row++) {
            col.push_back(repr);
        }
        col.push_back(fmt::format(notation, "..."));
        for (int row = size - rows / 2; row < size; row++) {
            col.push_back(repr);
        }
    }
    return col;
}

} // end unnamed namespace

string SolutionArray::info(const vector<string>& keys, int rows, int width)
{
    fmt::memory_buffer b;
    int col_width = 12; // targeted maximum column width
    vector<string> components;
    if (keys.size()) {
        components = keys;
    } else {
        components = componentNames();
    }
    try {
        // build columns
        vector<long int> index;
        for (const auto ix : m_active) {
            index.push_back(ix);
        }
        vector<vector<string>> cols = {integerColumn("", index, rows, col_width)};
        vector<vector<string>> tail; // trailing columns in reverse order
        size_t size = cols.back().size();

        // assemble columns fitting within a maximum width; if this width is exceeded,
        // a "..." separator is inserted close to the center. Accordingly, the matrix
        // needs to be assembled from two halves.
        int front = 0;
        int back = len(components) - 1;
        int fLen = len(cols.back()[0]);
        int bLen = 0;
        int sep = 5; // separator width
        bool done = false;
        while (!done && front <= back) {
            string key;
            while (bLen + sep <= fLen && front <= back) {
                // add trailing columns
                key = components[back];
                auto col = formatColumn(key, getComponent(key), rows, col_width);
                if (len(col[0]) + fLen + bLen + sep > width) {
                    done = true;
                    break;
                }
                tail.push_back(col);
                bLen += len(tail.back()[0]);
                back--;
            }
            if (done || front > back) {
                break;
            }
            while (fLen <= bLen + sep && front <= back) {
                // add leading columns
                key = components[front];
                auto col = formatColumn(key, getComponent(key), rows, col_width);
                if (len(col[0]) + fLen + bLen + sep > width) {
                    done = true;
                    break;
                }
                cols.push_back(col);
                fLen += len(cols.back()[0]);
                front++;
            }
        }
        if (cols.size() + tail.size() < components.size() + 1) {
            // add separator
            cols.push_back(vector<string>(size + 1, "  ..."));
        }
        // copy trailing columns
        cols.insert(cols.end(), tail.rbegin(), tail.rend());

        // assemble formatted output
        for (size_t row = 0; row < size; row++) {
            for (const auto& col : cols) {
                fmt_append(b, col[row]);
            }
            fmt_append(b, "\n");
        }

        // add size information
        fmt_append(b, "\n[{} rows x {} components; state='{}']",
            m_size, components.size(), m_sol->thermo()->nativeMode());

    } catch (CanteraError& err) {
        return to_string(b) + err.what();
    }
    return to_string(b);
}

shared_ptr<ThermoPhase> SolutionArray::thermo()
{
    return m_sol->thermo();
}

vector<string> SolutionArray::componentNames() const
{
    vector<string> components;
    // leading auxiliary components
    int pos = 0;
    while (m_order->count(pos)) {
        components.push_back(m_order->at(pos));
        pos++;
    }

    // state information
    auto phase = m_sol->thermo();
    for (auto code : phase->nativeMode()) {
        string name = string(1, code);
        if (name == "X" || name == "Y") {
            for (auto& spc : phase->speciesNames()) {
                components.push_back(spc);
            }
        } else {
            components.push_back(name);
        }
    }

    // trailing auxiliary components
    pos = -1;
    while (m_order->count(pos)) {
        components.push_back(m_order->at(pos));
        pos--;
    }

    return components;
}

void SolutionArray::addExtra(const string& name, bool back)
{
    if (m_extra->count(name)) {
        throw CanteraError("SolutionArray::addExtra",
            "Component '{}' already exists.", name);
    }
    (*m_extra)[name] = AnyValue();
    if (back) {
        if (m_order->size()) {
            // add name at end of back components
            m_order->emplace(m_order->begin()->first - 1, name);
        } else {
            // first name after state components
            m_order->emplace(-1, name);
        }
    } else {
        if (m_order->size()) {
            // insert name at end of front components
            m_order->emplace(m_order->rbegin()->first + 1, name);
        } else {
            // name in leading position
            m_order->emplace(0, name);
        }
    }
}

vector<string> SolutionArray::listExtra(bool all) const
{
    vector<string> names;
    int pos = 0;
    while (m_order->count(pos)) {
        const auto& name = m_order->at(pos);
        if (all || !m_extra->at(name).is<void>()) {
            names.push_back(name);
        }
        pos++;
    }

    // trailing auxiliary components
    pos = -1;
    while (m_order->count(pos)) {
        const auto& name = m_order->at(pos);
        if (all || !m_extra->at(name).is<void>()) {
            names.push_back(name);
        }
        pos--;
    }
    return names;
}

bool SolutionArray::hasComponent(const string& name) const
{
    if (m_extra->count(name)) {
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

AnyValue SolutionArray::getComponent(const string& name) const
{
    if (!hasComponent(name)) {
        throw CanteraError("SolutionArray::getComponent",
            "Unknown component '{}'.", name);
    }

    AnyValue out;
    if (m_extra->count(name)) {
        // extra component
        const auto& extra = m_extra->at(name);
        if (extra.is<void>()) {
            return AnyValue();
        }
        if (m_size == m_dataSize) {
            return extra; // slicing not necessary
        }
        if (extra.isVector<long int>()) {
            return getSingle<long int>(extra, m_active);
        }
        if (extra.isVector<double>()) {
            return getSingle<double>(extra, m_active);
        }
        if (extra.isVector<string>()) {
            return getSingle<string>(extra, m_active);
        }
        if (extra.isVector<vector<double>>()) {
            return getMulti<double>(extra, m_active);
        }
        if (extra.isVector<vector<long int>>()) {
            return getMulti<long int>(extra, m_active);
        }
        if (extra.isVector<vector<string>>()) {
            return getMulti<string>(extra, m_active);
        }
        throw NotImplementedError("SolutionArray::getComponent",
            "Unable to get sliced data for component '{}' with type '{}'.",
            name, extra.type_str());
    }

    // component is part of state information
    vector<double> data(m_size);
    size_t ix = m_sol->thermo()->speciesIndex(name);
    if (ix == npos) {
        // state other than species
        ix = m_sol->thermo()->nativeState()[name];
    } else {
        // species information
        ix += m_stride - m_sol->thermo()->nSpecies();
    }
    for (size_t k = 0; k < m_size; ++k) {
        data[k] = (*m_data)[m_active[k] * m_stride + ix];
    }
    out = data;
    return out;
}

bool isSimpleVector(const AnyValue& any) {
    return any.isVector<double>() || any.isVector<long int>() ||
        any.isVector<string>() || any.isVector<bool>() ||
        any.isVector<vector<double>>() || any.isVector<vector<long int>>() ||
        any.isVector<vector<string>>() || any.isVector<vector<bool>>();
}

void SolutionArray::setComponent(const string& name, const AnyValue& data)
{
    if (!hasComponent(name)) {
        throw CanteraError("SolutionArray::setComponent",
            "Unknown component '{}'.", name);
    }
    if (m_extra->count(name)) {
        _setExtra(name, data);
        return;
    }

    size_t size = data.vectorSize();
    if (size == npos) {
        throw CanteraError("SolutionArray::setComponent",
            "Invalid type of component '{}': expected simple array type, "
            "but received '{}'.", name, data.type_str());
    }
    if (size != m_size) {
        throw CanteraError("SolutionArray::setComponent",
            "Invalid size of component '{}': expected size {} but received {}.",
            name, m_size, size);
    }

    auto& vec = data.asVector<double>();
    size_t ix = m_sol->thermo()->speciesIndex(name);
    if (ix == npos) {
        ix = m_sol->thermo()->nativeState()[name];
    } else {
        ix += m_stride - m_sol->thermo()->nSpecies();
    }
    for (size_t k = 0; k < m_size; ++k) {
        (*m_data)[m_active[k] * m_stride + ix] = vec[k];
    }
}

void SolutionArray::setLoc(int loc, bool restore)
{
    size_t loc_ = static_cast<size_t>(loc);
    if (m_size == 0) {
        throw CanteraError("SolutionArray::setLoc",
            "Unable to set location in empty SolutionArray.");
    } else if (loc < 0) {
        if (m_loc == npos) {
            throw CanteraError("SolutionArray::setLoc",
                "Both current and buffered indices are invalid.");
        }
        return;
    } else if (static_cast<size_t>(m_active[loc_]) == m_loc) {
        return;
    } else if (loc_ >= m_size) {
        throw IndexError("SolutionArray::setLoc", "indices", loc_, m_size - 1);
    }
    m_loc = static_cast<size_t>(m_active[loc_]);
    if (restore) {
        size_t nState = m_sol->thermo()->stateSize();
        m_sol->thermo()->restoreState(nState, m_data->data() + m_loc * m_stride);
    }
}

void SolutionArray::updateState(int loc)
{
    setLoc(loc, false);
    size_t nState = m_sol->thermo()->stateSize();
    m_sol->thermo()->saveState(nState, m_data->data() + m_loc * m_stride);
}

vector<double> SolutionArray::getState(int loc)
{
    setLoc(loc);
    size_t nState = m_sol->thermo()->stateSize();
    vector<double> out(nState);
    m_sol->thermo()->saveState(out); // thermo contains current state
    return out;
}

void SolutionArray::setState(int loc, const vector<double>& state)
{
    size_t nState = m_sol->thermo()->stateSize();
    if (state.size() != nState) {
        throw CanteraError("SolutionArray::setState",
            "Expected array to have length {}, but received an array of length {}.",
            nState, state.size());
    }
    setLoc(loc, false);
    m_sol->thermo()->restoreState(state);
    m_sol->thermo()->saveState(nState, m_data->data() + m_loc * m_stride);
}

void SolutionArray::normalize() {
    auto phase = m_sol->thermo();
    auto nativeState = phase->nativeState();
    if (nativeState.size() < 3) {
        return;
    }
    size_t nState = phase->stateSize();
    vector<double> out(nState);
    if (nativeState.count("Y")) {
        size_t offset = nativeState["Y"];
        for (int loc = 0; loc < static_cast<int>(m_size); loc++) {
            setLoc(loc, true); // set location and restore state
            phase->setMassFractions(m_data->data() + m_loc * m_stride + offset);
            m_sol->thermo()->saveState(out);
            setState(loc, out);
        }
    } else if (nativeState.count("X")) {
        size_t offset = nativeState["X"];
        for (int loc = 0; loc < static_cast<int>(m_size); loc++) {
            setLoc(loc, true); // set location and restore state
            phase->setMoleFractions(m_data->data() + m_loc * m_stride + offset);
            m_sol->thermo()->saveState(out);
            setState(loc, out);
        }
    } else {
        throw NotImplementedError("SolutionArray::normalize",
            "Not implemented for mode '{}'.", phase->nativeMode());
    }
}

AnyMap SolutionArray::getAuxiliary(int loc)
{
    setLoc(loc);
    AnyMap out;
    for (const auto& [key, extra] : *m_extra) {
        if (extra.is<void>()) {
            out[key] = extra;
        } else if (extra.isVector<long int>()) {
            out[key] = extra.asVector<long int>()[m_loc];
        } else if (extra.isVector<double>()) {
            out[key] = extra.asVector<double>()[m_loc];
        } else if (extra.isVector<string>()) {
            out[key] = extra.asVector<string>()[m_loc];
        } else if (extra.isVector<vector<long int>>()) {
            out[key] = extra.asVector<vector<long int>>()[m_loc];
        } else if (extra.isVector<vector<double>>()) {
            out[key] = extra.asVector<vector<double>>()[m_loc];
        } else if (extra.isVector<vector<string>>()) {
            out[key] = extra.asVector<vector<string>>()[m_loc];
        } else {
            throw NotImplementedError("SolutionArray::getAuxiliary",
                "Unable to retrieve data for component '{}' with type '{}'.",
                key, extra.type_str());
        }
    }
    return out;
}

void SolutionArray::setAuxiliary(int loc, const AnyMap& data)
{
    setLoc(loc, false);
    for (const auto& [name, value] : data) {
        if (!m_extra->count(name)) {
            throw CanteraError("SolutionArray::setAuxiliary",
                "Unknown auxiliary component '{}'.", name);
        }
        auto& extra = m_extra->at(name);
        if (extra.is<void>()) {
            if (m_dataSize > 1) {
                throw CanteraError("SolutionArray::setAuxiliary",
                    "Unable to set location for type '{}': "
                    "component is not initialized.", name);
            }
            _initExtra(name, value);
            _resizeExtra(name);
        }
        try {
            if (extra.isVector<long int>()) {
                setAuxiliarySingle<long int>(m_loc, extra, value);
            } else if (extra.isVector<double>()) {
                setAuxiliarySingle<double>(m_loc, extra, value);
            } else if (extra.isVector<string>()) {
                setAuxiliarySingle<string>(m_loc, extra, value);
            } else if (extra.isVector<vector<long int>>()) {
                setAuxiliaryMulti<long int>(m_loc, extra, value);
            } else if (extra.isVector<vector<double>>()) {
                setAuxiliaryMulti<double>(m_loc, extra, value);
            } else if (extra.isVector<vector<string>>()) {
                setAuxiliaryMulti<string>(m_loc, extra, value);
            } else {
                throw CanteraError("SolutionArray::setAuxiliary",
                    "Unable to set entry for type '{}'.", extra.type_str());
            }
        } catch (CanteraError& err) {
            // make failed type conversions traceable
            throw CanteraError("SolutionArray::setAuxiliary",
                "Encountered incompatible value for component '{}':\n{}",
                name, err.getMessage());
        }
    }
}

AnyMap preamble(const string& desc)
{
    AnyMap data;
    if (desc.size()) {
        data["description"] = desc;
    }
    data["generator"] = "Cantera SolutionArray";
    data["cantera-version"] = CANTERA_VERSION;
    // escape commit to ensure commits are read correctly from YAML
    // example: prevent '3491027e7' from being converted to an integer
    data["git-commit"] = "'" + gitCommit() + "'";

    // Add a timestamp indicating the current time
    time_t aclock;
    ::time(&aclock); // Get time in seconds
    struct tm* newtime = localtime(&aclock); // Convert time to struct tm form
    data["date"] = stripnonprint(asctime(newtime));

    // Force metadata fields to the top of the file
    if (data.hasKey("description")) {
        data["description"].setLoc(-6, 0);
    }
    data["generator"].setLoc(-5, 0);
    data["cantera-version"].setLoc(-4, 0);
    data["git-commit"].setLoc(-3, 0);
    data["date"].setLoc(-2, 0);

    return data;
}

AnyMap& openField(AnyMap& root, const string& name)
{
    if (!name.size()) {
        return root;
    }

    // locate field based on 'name'
    vector<string> tokens;
    tokenizePath(name, tokens);
    AnyMap* ptr = &root; // use raw pointer to avoid copying
    string path = "";
    for (auto& field : tokens) {
        path += "/" + field;
        AnyMap& sub = *ptr;
        if (sub.hasKey(field) && !sub[field].is<AnyMap>()) {
            throw CanteraError("openField",
                "Encountered invalid existing field '{}'.", path);
        } else if (!sub.hasKey(field)) {
            sub[field] = AnyMap();
        }
        ptr = &sub[field].as<AnyMap>();
    }
    return *ptr;
}

void SolutionArray::writeHeader(const string& fname, const string& name,
                                const string& desc, bool overwrite)
{
    Storage file(fname, true);
    if (file.checkGroup(name, true)) {
        if (!overwrite) {
            throw CanteraError("SolutionArray::writeHeader",
                "Group name '{}' exists; use 'overwrite' argument to overwrite.", name);
        }
        file.deleteGroup(name);
        file.checkGroup(name, true);
    }
    file.writeAttributes(name, preamble(desc));
}

void SolutionArray::writeHeader(AnyMap& root, const string& name,
                                const string& desc, bool overwrite)
{
    AnyMap& data = openField(root, name);
    if (!data.empty() && !overwrite) {
        throw CanteraError("SolutionArray::writeHeader",
            "Field name '{}' exists; use 'overwrite' argument to overwrite.", name);
    }
    data.update(preamble(desc));
}

void SolutionArray::writeEntry(const string& fname, bool overwrite, const string& basis)
{
    if (apiNdim() != 1) {
        throw CanteraError("SolutionArray::writeEntry",
            "Tabular output of CSV data only works for 1D SolutionArray objects.");
    }
    set<string> speciesNames;
    for (const auto& species : m_sol->thermo()->speciesNames()) {
        speciesNames.insert(species);
    }
    bool mole;
    if (basis == "") {
        const auto& nativeState = m_sol->thermo()->nativeState();
        mole = nativeState.find("X") != nativeState.end();
    } else if (basis == "X" || basis == "mole") {
        mole = true;
    } else if (basis == "Y" || basis == "mass") {
        mole = false;
    } else {
        throw CanteraError("SolutionArray::writeEntry",
            "Invalid species basis '{}'.", basis);
    }

    auto names = componentNames();
    size_t last = names.size() - 1;
    vector<AnyValue> components;
    vector<bool> isSpecies;
    fmt::memory_buffer header;
    for (const auto& key : names) {
        string label = key;
        size_t col;
        if (speciesNames.find(key) == speciesNames.end()) {
            // Pre-read component vectors
            isSpecies.push_back(false);
            components.emplace_back(getComponent(key));
            col = components.size() - 1;
            if (!components[col].isVector<double>() &&
                !components[col].isVector<long int>() &&
                !components[col].isVector<string>())
            {
                throw CanteraError("SolutionArray::writeEntry",
                    "Multi-dimensional column '{}' is not supported for CSV output.",
                    key);
            }
        } else {
            // Delay reading species data as base can be either mole or mass
            isSpecies.push_back(true);
            components.emplace_back(AnyValue());
            col = components.size() - 1;
            if (mole) {
                label = "X_" + label;
            } else {
                label = "Y_" + label;
            }
        }
        if (label.find("\"") != string::npos || label.find("\n") != string::npos) {
            throw NotImplementedError("SolutionArray::writeEntry",
                "Detected column name containing double quotes or line feeds: '{}'.",
                label);
        }
        string sep = (col == last) ? "" : ",";
        if (label.find(",") != string::npos) {
            fmt_append(header, "\"{}\"{}", label, sep);
        } else {
            fmt_append(header, "{}{}", label, sep);
        }
    }

    // (Most) potential exceptions have been thrown; start writing data to file
    if (std::ifstream(fname).good()) {
        if (!overwrite) {
            throw CanteraError("SolutionArray::writeEntry",
                "File '{}' already exists; use option 'overwrite' to replace CSV file.",
                fname);
        }
        std::remove(fname.c_str());
    }
    std::ofstream output(fname);
    output << to_string(header) << std::endl;

    size_t maxLen = npos;
    vector<double> buf(speciesNames.size(), 0.);
    for (int row = 0; row < static_cast<int>(m_size); row++) {
        fmt::memory_buffer line;
        if (maxLen != npos) {
            line.reserve(maxLen);
        }
        setLoc(row);
        if (mole) {
            m_sol->thermo()->getMoleFractions(buf.data());
        } else {
            m_sol->thermo()->getMassFractions(buf.data());
        }

        size_t idx = 0;
        for (size_t col = 0; col < components.size(); col++) {
            string sep = (col == last) ? "" : ",";
            if (isSpecies[col]) {
                fmt_append(line, "{:.9g}{}", buf[idx++], sep);
            } else {
                auto& data = components[col];
                if (data.isVector<double>()) {
                    fmt_append(line, "{:.9g}{}", data.asVector<double>()[row], sep);
                } else if (data.isVector<long int>()) {
                    fmt_append(line, "{}{}", data.asVector<long int>()[row], sep);
                } else if (data.isVector<bool>()) {
                    fmt_append(line, "{}{}",
                               static_cast<bool>(data.asVector<bool>()[row]), sep);
                } else {
                    auto value = data.asVector<string>()[row];
                    if (value.find("\"") != string::npos ||
                        value.find("\n") != string::npos)
                    {
                        throw NotImplementedError("SolutionArray::writeEntry",
                            "Detected value containing double quotes or line feeds: "
                            "'{}'", value);
                    }
                    if (value.find(",") != string::npos) {
                        fmt_append(line, "\"{}\"{}", value, sep);
                    } else {
                        fmt_append(line, "{}{}", value, sep);
                    }
                }
            }
        }
        output << to_string(line) << std::endl;
        maxLen = std::max(maxLen, line.size());
    }
}

void SolutionArray::writeEntry(const string& fname, const string& name,
                               const string& sub, bool overwrite, int compression)
{
    if (name == "") {
        throw CanteraError("SolutionArray::writeEntry",
            "Group name specifying root location must not be empty.");
    }
    if (m_size < m_dataSize) {
        throw NotImplementedError("SolutionArray::writeEntry",
            "Unable to save sliced data.");
    }
    Storage file(fname, true);
    if (compression) {
        file.setCompressionLevel(compression);
    }
    string path = name;
    if (sub != "") {
        path += "/" + sub;
    } else {
        path += "/data";
    }
    if (file.checkGroup(path, true)) {
        if (!overwrite) {
            throw CanteraError("SolutionArray::writeEntry",
                "Group name '{}' exists; use 'overwrite' argument to overwrite.", name);
        }
        file.deleteGroup(path);
        file.checkGroup(path, true);
    }
    file.writeAttributes(path, m_meta);
    AnyMap more;
    if (apiNdim() == 1) {
        more["size"] = int(m_dataSize);
    } else {
        more["api-shape"] = m_apiShape;
    }
    more["components"] = componentNames();
    file.writeAttributes(path, more);
    if (!m_dataSize) {
        return;
    }

    const auto& nativeState = m_sol->thermo()->nativeState();
    size_t nSpecies = m_sol->thermo()->nSpecies();
    for (auto& [key, offset] : nativeState) {
        if (key == "X" || key == "Y") {
            vector<vector<double>> prop;
            for (size_t i = 0; i < m_size; i++) {
                size_t first = offset + i * m_stride;
                prop.emplace_back(m_data->begin() + first,
                                  m_data->begin() + first + nSpecies);
            }
            AnyValue data;
            data = prop;
            file.writeData(path, key, data);
        } else {
            auto data = getComponent(key);
            file.writeData(path, key, data);
        }
    }

    for (const auto& [key, value] : *m_extra) {
        if (isSimpleVector(value)) {
            file.writeData(path, key, value);
        } else if (value.is<void>()) {
            // skip unintialized component
        } else {
            throw NotImplementedError("SolutionArray::writeEntry",
                "Unable to save component '{}' with data type {}.",
                key, value.type_str());
        }
    }
}

void SolutionArray::writeEntry(AnyMap& root, const string& name, const string& sub,
                               bool overwrite)
{
    if (name == "") {
        throw CanteraError("SolutionArray::writeEntry",
            "Field name specifying root location must not be empty.");
    }
    if (m_size < m_dataSize) {
        throw NotImplementedError("SolutionArray::writeEntry",
            "Unable to save sliced data.");
    }
    string path = name;
    if (sub != "") {
        path += "/" + sub;
    } else {
        path += "/data";
    }
    AnyMap& data = openField(root, path);
    bool preexisting = !data.empty();
    if (preexisting && !overwrite) {
        throw CanteraError("SolutionArray::writeEntry",
            "Field name '{}' exists; use 'overwrite' argument to overwrite.", name);
    }
    if (apiNdim() == 1) {
        data["size"] = int(m_dataSize);
    } else {
        data["api-shape"] = m_apiShape;
    }
    data.update(m_meta);

    for (auto& [key, value] : *m_extra) {
        data[key] = value;
    }

    auto phase = m_sol->thermo();
    if (m_size == 1) {
        setLoc(0);
        data["temperature"] = phase->temperature();
        data["pressure"] = phase->pressure();
        auto surf = std::dynamic_pointer_cast<SurfPhase>(phase);
        auto nSpecies = phase->nSpecies();
        vector<double> values(nSpecies);
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
        for (auto& [key, offset] : nativeState) {
            if (key == "X" || key == "Y") {
                for (auto& spc : phase->speciesNames()) {
                    data[spc] = getComponent(spc);
                }
                data["basis"] = key == "X" ? "mole" : "mass";
            } else {
                data[key] = getComponent(key);
            }
        }
        data["components"] = componentNames();
    }

    // If this is not replacing an existing solution, put it at the end
    if (!preexisting) {
        data.setLoc(INT_MAX, 0);
    }
}

void SolutionArray::append(const vector<double>& state, const AnyMap& extra)
{
    if (apiNdim() > 1) {
        throw NotImplementedError("SolutionArray::append",
            "Unable to append multi-dimensional arrays.");
    }

    int pos = size();
    resize(pos + 1);
    try {
        setState(pos, state);
        setAuxiliary(pos, extra);
    } catch (CanteraError& err) {
        // undo resize and rethrow error
        resize(pos);
        throw CanteraError("SolutionArray::append", err.getMessage());
    }
}

void SolutionArray::save(const string& fname, const string& name, const string& sub,
                         const string& desc, bool overwrite, int compression,
                         const string& basis)
{
    if (m_size < m_dataSize) {
        throw NotImplementedError("SolutionArray::save",
                                  "Unable to save sliced data.");
    }
    size_t dot = fname.find_last_of(".");
    string extension = (dot != npos) ? toLowerCopy(fname.substr(dot + 1)) : "";
    if (extension == "csv") {
        if (name != "") {
            warn_user("SolutionArray::save",
                      "Parameter 'name' not used for CSV output.");
        }
        writeEntry(fname, overwrite, basis);
        return;
    }
    if (basis != "") {
        warn_user("SolutionArray::save",
            "Argument 'basis' is not used for HDF or YAML output.", basis);
    }
    if (extension == "h5" || extension == "hdf"  || extension == "hdf5") {
        writeHeader(fname, name, desc, overwrite);
        writeEntry(fname, name, sub, true, compression);
        return;
    }
    if (extension == "yaml" || extension == "yml") {
        // Check for an existing file and load it if present
        AnyMap data;
        if (std::ifstream(fname).good()) {
            data = AnyMap::fromYamlFile(fname);
        }
        writeHeader(data, name, desc, overwrite);
        writeEntry(data, name, sub, true);

        // Write the output file and remove the now-outdated cached file
        std::ofstream out(fname);
        out << data.toYamlString();
        AnyMap::clearCachedFile(fname);
        return;
    }
    throw CanteraError("SolutionArray::save",
                       "Unknown file extension '{}'.", extension);
}

AnyMap SolutionArray::readHeader(const string& fname, const string& name)
{
    Storage file(fname, false);
    file.checkGroup(name);
    return file.readAttributes(name, false);
}

const AnyMap& locateField(const AnyMap& root, const string& name)
{
    if (!name.size()) {
        return root;
    }

    // locate field based on 'name'
    vector<string> tokens;
    tokenizePath(name, tokens);
    const AnyMap* ptr = &root; // use raw pointer to avoid copying
    string path = "";
    for (auto& field : tokens) {
        path += "/" + field;
        const AnyMap& sub = *ptr;
        if (!sub.hasKey(field) || !sub[field].is<AnyMap>()) {
            throw CanteraError("SolutionArray::locateField",
                "No field or solution with name '{}'.", path);
        }
        ptr = &sub[field].as<AnyMap>();
    }
    return *ptr;
}

AnyMap SolutionArray::readHeader(const AnyMap& root, const string& name)
{
    auto sub = locateField(root, name);
    AnyMap header;
    for (const auto& [key, value] : sub) {
        if (!sub[key].is<AnyMap>()) {
            header[key] = value;
        }
    }
    return header;
}

AnyMap SolutionArray::restore(const string& fname,
                              const string& name, const string& sub)
{
    size_t dot = fname.find_last_of(".");
    string extension = (dot != npos) ? toLowerCopy(fname.substr(dot + 1)) : "";
    AnyMap header;
    if (extension == "csv") {
        throw NotImplementedError("SolutionArray::restore",
            "CSV import not implemented; if using Python, data can be imported via "
            "'read_csv' instead.");
    }
    if (extension == "h5" || extension == "hdf"  || extension == "hdf5") {
        readEntry(fname, name, sub);
        header = readHeader(fname, name);
    } else if (extension == "yaml" || extension == "yml") {
        const AnyMap& root = AnyMap::fromYamlFile(fname);
        readEntry(root, name, sub);
        header = readHeader(root, name);
    } else {
        throw CanteraError("SolutionArray::restore",
            "Unknown file extension '{}'; supported extensions include "
            "'h5'/'hdf'/'hdf5' and 'yml'/'yaml'.", extension);
    }
    return header;
}

void SolutionArray::_initExtra(const string& name, const AnyValue& value)
{
    if (!m_extra->count(name)) {
        throw CanteraError("SolutionArray::_initExtra",
            "Component '{}' does not exist.", name);
    }
    auto& extra = (*m_extra)[name];
    if (!extra.is<void>()) {
        throw CanteraError("SolutionArray::_initExtra",
            "Component '{}' is already initialized.", name);
    }
    try {
        if (value.is<long int>()) {
            extra = vector<long int>(m_dataSize, value.as<long int>());
        } else if (value.is<double>()) {
            extra = vector<double>(m_dataSize, value.as<double>());
        } else if (value.is<string>()) {
            extra = vector<string>(m_dataSize, value.as<string>());
        } else if (value.isVector<long int>()) {
            extra = vector<vector<long int>>(m_dataSize, value.asVector<long int>());
        } else if (value.isVector<double>()) {
            extra = vector<vector<double>>(m_dataSize, value.asVector<double>());
        } else if (value.isVector<string>()) {
            extra = vector<vector<string>>(m_dataSize, value.asVector<string>());
        } else if (value.is<void>()) {
            // placeholder for unknown type; settable in setComponent
            extra = AnyValue();
        } else {
            throw NotImplementedError("SolutionArray::_initExtra",
                "Unable to initialize component '{}' with type '{}'.",
                name, value.type_str());
        }
    } catch (CanteraError& err) {
        // make failed type conversions traceable
        throw CanteraError("SolutionArray::_initExtra",
            "Encountered incompatible value for initializing component '{}':\n{}",
            name, err.getMessage());
    }
}

void SolutionArray::_resizeExtra(const string& name, const AnyValue& value)
{
    if (!m_extra->count(name)) {
        throw CanteraError("SolutionArray::_resizeExtra",
            "Component '{}' does not exist.", name);
    }
    auto& extra = (*m_extra)[name];
    if (extra.is<void>()) {
        // cannot resize placeholder (uninitialized component)
        return;
    }

    try {
        if (extra.isVector<long int>()) {
            resizeSingle<long int>(extra, m_dataSize, value);
        } else if (extra.isVector<double>()) {
            resizeSingle<double>(extra, m_dataSize, value);
        } else if (extra.isVector<string>()) {
            resizeSingle<string>(extra, m_dataSize, value);
        } else if (extra.isVector<vector<double>>()) {
            resizeMulti<double>(extra, m_dataSize, value);
        } else if (extra.isVector<vector<long int>>()) {
            resizeMulti<long int>(extra, m_dataSize, value);
        } else if (extra.isVector<vector<string>>()) {
            resizeMulti<string>(extra, m_dataSize, value);
        } else {
            throw NotImplementedError("SolutionArray::_resizeExtra",
                "Unable to resize using type '{}'.", extra.type_str());
        }
    } catch (CanteraError& err) {
        // make failed type conversions traceable
        throw CanteraError("SolutionArray::_resizeExtra",
            "Encountered incompatible value for resizing component '{}':\n{}",
            name, err.getMessage());
    }
}

void SolutionArray::_setExtra(const string& name, const AnyValue& data)
{
    if (!m_extra->count(name)) {
        throw CanteraError("SolutionArray::_setExtra",
                           "Extra component '{}' does not exist.", name);
    }

    auto& extra = m_extra->at(name);
    if (data.is<void>() && m_size == m_dataSize) {
        // reset placeholder
        extra = AnyValue();
        return;
    }

    // uninitialized component
    if (extra.is<void>()) {
        if (m_size != m_dataSize) {
            throw CanteraError("SolutionArray::_setExtra",
                "Unable to replace '{}' for sliced data.", name);
        }
        if (data.vectorSize() == m_dataSize || data.matrixShape().first == m_dataSize) {
            // assign entire component
            extra = data;
            return;
        }
        // data contains default value
        if (m_dataSize == 0 && (data.isScalar() || data.vectorSize() > 0)) {
            throw CanteraError("SolutionArray::_setExtra",
                "Unable to initialize '{}' with non-empty values when SolutionArray is "
                "empty.", name);
        }
        if (m_dataSize && data.vectorSize() == 0) {
            throw CanteraError("SolutionArray::_setExtra",
                "Unable to initialize '{}' with empty array when SolutionArray is not "
                "empty." , name);
        }
        _initExtra(name, data);
        _resizeExtra(name, data);
        return;
    }

    if (data.is<long int>()) {
        setScalar<long int>(extra, data, m_active);
    } else if (data.is<double>()) {
        setScalar<double>(extra, data, m_active);
    } else if (data.is<string>()) {
        setScalar<string>(extra, data, m_active);
    } else if (data.isVector<long int>()) {
        setSingle<long int>(extra, data, m_active);
    } else if (data.isVector<double>()) {
        setSingle<double>(extra, data, m_active);
    } else if (data.isVector<string>()) {
        setSingle<string>(extra, data, m_active);
    } else if (data.isVector<vector<long int>>()) {
        setMulti<long int>(extra, data, m_active);
    } else if (data.isVector<vector<double>>()) {
        setMulti<double>(extra, data, m_active);
    } else if (data.isVector<vector<string>>()) {
        setMulti<string>(extra, data, m_active);
    } else {
        throw NotImplementedError("SolutionArray::_setExtra",
            "Unable to set sliced data for component '{}' with type '{}'.",
            name, data.type_str());
    }
}

string SolutionArray::_detectMode(const set<string>& names, bool native)
{
    // check set of available names against state acronyms defined by Phase::fullStates
    string mode = "";
    const auto& nativeState = m_sol->thermo()->nativeState();
    bool usesNativeState = false;
    auto surf = std::dynamic_pointer_cast<SurfPhase>(m_sol->thermo());
    bool found;
    for (const auto& item : m_sol->thermo()->fullStates()) {
        found = true;
        string name;
        usesNativeState = true;
        for (size_t i = 0; i < item.size(); i++) {
            // pick i-th letter from "full" state acronym
            name = string(1, item[i]);
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
    if (!found) {
        if (surf && names.count("T") && names.count("X") && names.count("density")) {
            // Legacy HDF format erroneously uses density (Cantera < 3.0)
            return "legacySurf";
        }
        if (names.count("mass-flux") && names.count("mass-fractions")) {
            // Legacy YAML format stores incomplete state information (Cantera < 3.0)
            return "legacyInlet";
        }
        throw CanteraError("SolutionArray::_detectMode",
            "Detected incomplete thermodynamic information. Full states for a '{}' "
            "phase\nmay be defined by the following modes:\n'{}'\n"
            "Available data are: '{}'", m_sol->thermo()->type(),
            ba::join(m_sol->thermo()->fullStates(), "', '"), ba::join(names, "', '"));
    }
    if (usesNativeState && native) {
        return "native";
    }
    return mode;
}

set<string> SolutionArray::_stateProperties(
    const string& mode, bool alias)
{
    set<string> states;
    if (mode == "native") {
        for (const auto& [name, offset] : m_sol->thermo()->nativeState()) {
            states.insert(alias ? aliasMap.at(name) : name);
        }
    } else {
        for (const auto& m : mode) {
            const string name = string(1, m);
            states.insert(alias ? aliasMap.at(name) : name);
        }
    }

    return states;
}

string getName(const set<string>& names, const string& name)
{
    if (names.count(name)) {
        return name;
    }
    const auto& alias = aliasMap.at(name);
    if (names.count(alias)) {
        return alias;
    }
    return name; // let exception be thrown elsewhere
}

void SolutionArray::readEntry(const string& fname, const string& name,
                              const string& sub)
{
    Storage file(fname, false);
    if (name == "") {
        throw CanteraError("SolutionArray::readEntry",
            "Group name specifying root location must not be empty.");
    }
    string path = name;
    if (sub != "" && file.checkGroup(name + "/" + sub, true)) {
        path += "/" + sub;
    } else if (sub == "" && file.checkGroup(name + "/data", true)) {
        // default data location
        path += "/data";
    }
    if (!file.checkGroup(path)) {
        throw CanteraError("SolutionArray::readEntry",
            "Group name specifying data entry is empty.");
    }
    m_extra->clear();
    auto [size, names] = file.contents(path);
    m_meta = file.readAttributes(path, true);
    if (m_meta.hasKey("size")) {
        // one-dimensional array
        resize(m_meta["size"].as<long int>());
        m_meta.erase("size");
    } else if (m_meta.hasKey("api-shape")) {
        // API uses multiple dimensions to interpret C++ SolutionArray
        setApiShape(m_meta["api-shape"].asVector<long int>());
        m_meta.erase("api-shape");
    } else {
        // legacy format; size needs to be detected
        resize(static_cast<int>(size));
    }

    if (m_size == 0) {
        return;
    }

    // determine storage mode of state data
    string mode = _detectMode(names);
    set<string> states = _stateProperties(mode);
    if (states.count("C")) {
        if (names.count("X")) {
            states.erase("C");
            states.insert("X");
            mode = mode.substr(0, 2) +  "X";
        } else if (names.count("Y")) {
            states.erase("C");
            states.insert("Y");
            mode = mode.substr(0, 2) +  "Y";
        }
    }

    // restore state data
    size_t nSpecies = m_sol->thermo()->nSpecies();
    size_t nState = m_sol->thermo()->stateSize();
    const auto& nativeStates = m_sol->thermo()->nativeState();
    if (mode == "native") {
        // native state can be written directly into data storage
        for (const auto& [name, offset] : nativeStates) {
            if (name == "X" || name == "Y") {
                AnyValue data;
                data = file.readData(path, name, m_size, nSpecies);
                auto prop = data.asVector<vector<double>>();
                for (size_t i = 0; i < m_dataSize; i++) {
                    std::copy(prop[i].begin(), prop[i].end(),
                              m_data->data() + offset + i * m_stride);
                }
            } else {
                AnyValue data;
                data = file.readData(path, getName(names, name), m_dataSize, 0);
                setComponent(name, data);
            }
        }
    } else if (mode == "TPX") {
        AnyValue data;
        data = file.readData(path, getName(names, "T"), m_dataSize, 0);
        vector<double> T = std::move(data.asVector<double>());
        data = file.readData(path, getName(names, "P"), m_dataSize, 0);
        vector<double> P = std::move(data.asVector<double>());
        data = file.readData(path, "X", m_dataSize, nSpecies);
        vector<vector<double>> X = std::move(data.asVector<vector<double>>());
        for (size_t i = 0; i < m_dataSize; i++) {
            m_sol->thermo()->setMoleFractions_NoNorm(X[i].data());
            m_sol->thermo()->setState_TP(T[i], P[i]);
            m_sol->thermo()->saveState(nState, m_data->data() + i * m_stride);
        }
    } else if (mode == "TDX") {
        AnyValue data;
        data = file.readData(path, getName(names, "T"), m_dataSize, 0);
        vector<double> T = std::move(data.asVector<double>());
        data = file.readData(path, getName(names, "D"), m_dataSize, 0);
        vector<double> D = std::move(data.asVector<double>());
        data = file.readData(path, "X", m_dataSize, nSpecies);
        vector<vector<double>> X = std::move(data.asVector<vector<double>>());
        for (size_t i = 0; i < m_dataSize; i++) {
            m_sol->thermo()->setMoleFractions_NoNorm(X[i].data());
            m_sol->thermo()->setState_TD(T[i], D[i]);
            m_sol->thermo()->saveState(nState, m_data->data() + i * m_stride);
        }
    } else if (mode == "TPY") {
        AnyValue data;
        data = file.readData(path, getName(names, "T"), m_dataSize, 0);
        vector<double> T = std::move(data.asVector<double>());
        data = file.readData(path, getName(names, "P"), m_dataSize, 0);
        vector<double> P = std::move(data.asVector<double>());
        data = file.readData(path, "Y", m_dataSize, nSpecies);
        vector<vector<double>> Y = std::move(data.asVector<vector<double>>());
        for (size_t i = 0; i < m_dataSize; i++) {
            m_sol->thermo()->setMassFractions_NoNorm(Y[i].data());
            m_sol->thermo()->setState_TP(T[i], P[i]);
            m_sol->thermo()->saveState(nState, m_data->data() + i * m_stride);
        }
    } else if (mode == "legacySurf") {
        // erroneous TDX mode (should be TPX or TPY) - Sim1D (Cantera 2.5)
        AnyValue data;
        data = file.readData(path, getName(names, "T"), m_dataSize, 0);
        vector<double> T = std::move(data.asVector<double>());
        data = file.readData(path, "X", m_dataSize, nSpecies);
        vector<vector<double>> X = std::move(data.asVector<vector<double>>());
        for (size_t i = 0; i < m_dataSize; i++) {
            m_sol->thermo()->setMoleFractions_NoNorm(X[i].data());
            m_sol->thermo()->setTemperature(T[i]);
            m_sol->thermo()->saveState(nState, m_data->data() + i * m_stride);
        }
        warn_user("SolutionArray::readEntry",
            "Detected legacy HDF format with incomplete state information\nfor name "
            "'{}' (pressure missing).", path);
    } else if (mode == "") {
        throw CanteraError("SolutionArray::readEntry",
            "Data are not consistent with full state modes.");
    } else {
        throw NotImplementedError("SolutionArray::readEntry",
            "Import of '{}' data is not supported.", mode);
    }

    // restore remaining data
    if (m_meta.hasKey("components")) {
        const auto& components = m_meta["components"].asVector<string>();
        bool back = false;
        for (const auto& name : components) {
            if (hasComponent(name) || name == "X" || name == "Y") {
                back = true;
            } else {
                addExtra(name, back);
                AnyValue data;
                data = file.readData(path, name, m_dataSize);
                setComponent(name, data);
            }
        }
        m_meta.erase("components");
    } else {
        // data format used by Python h5py export (Cantera 2.5)
        warn_user("SolutionArray::readEntry", "Detected legacy HDF format.");
        for (const auto& name : names) {
            if (!hasComponent(name) && name != "X" && name != "Y") {
                addExtra(name);
                AnyValue data;
                data = file.readData(path, name, m_dataSize);
                setComponent(name, data);
            }
        }
    }
}

void SolutionArray::readEntry(const AnyMap& root, const string& name, const string& sub)
{
    if (name == "") {
        throw CanteraError("SolutionArray::readEntry",
            "Field name specifying root location must not be empty.");
    }
    auto path = locateField(root, name);
    if (path.hasKey("generator") && sub != "") {
        // read entry from subfolder (since Cantera 3.0)
        path = locateField(root, name + "/" + sub);
    } else if (sub == "" && path.hasKey("data")) {
        // default data location
        path = locateField(root, name + "/data");
    }

    // set size and initialize
    long size = 0;
    if (path.hasKey("size")) {
        // one-dimensional array
        resize(path["size"].asInt());
    } else if (path.hasKey("api-shape")) {
        // API uses multiple dimensions to interpret C++ SolutionArray
        auto& shape = path["api-shape"].asVector<long int>();
        setApiShape(shape);
    } else {
        // legacy format (Cantera 2.6)
        size = path.getInt("points", 0);
        if (!path.hasKey("T") && !path.hasKey("temperature")) {
            // overwrite size - Sim1D erroneously assigns '1'
            size = (long)0;
        }
        resize(static_cast<int>(size));
    }
    m_extra->clear();

    // restore data
    set<string> exclude = {"size", "api-shape", "points", "X", "Y"};
    set<string> names = path.keys();
    size_t nState = m_sol->thermo()->stateSize();
    if (m_dataSize == 0) {
        // no data points
    } else if (m_dataSize == 1) {
        // single data point
        string mode = _detectMode(names, false);
        if (mode == "TPY") {
            double T = path[getName(names, "T")].asDouble();
            double P = path[getName(names, "P")].asDouble();
            auto Y = path["mass-fractions"].asMap<double>();
            m_sol->thermo()->setState_TPY(T, P, Y);
        } else if (mode == "TPC") {
            auto surf = std::dynamic_pointer_cast<SurfPhase>(m_sol->thermo());
            double T = path[getName(names, "T")].asDouble();
            double P = path["pressure"].asDouble();
            m_sol->thermo()->setState_TP(T, P);
            auto cov = path["coverages"].asMap<double>();
            surf->setCoveragesByName(cov);
        } else if (mode == "legacyInlet") {
            // missing property - Sim1D (Cantera 2.6)
            mode = "TPY";
            double T = path[getName(names, "T")].asDouble();
            auto Y = path["mass-fractions"].asMap<double>();
            m_sol->thermo()->setState_TPY(T, m_sol->thermo()->pressure(), Y);
            warn_user("SolutionArray::readEntry",
                "Detected legacy YAML format with incomplete state information\n"
                "for name '{}' (pressure missing).", name + "/" + sub);
        } else if (mode == "") {
            throw CanteraError("SolutionArray::readEntry",
                "Data are not consistent with full state modes.");
        } else {
            throw NotImplementedError("SolutionArray::readEntry",
                "Import of '{}' data is not supported.", mode);
        }
        m_sol->thermo()->saveState(nState, m_data->data());
        auto props = _stateProperties(mode, true);
        exclude.insert(props.begin(), props.end());
    } else {
        // multiple data points
        if (path.hasKey("components")) {
            const auto& components = path["components"].asVector<string>();
            bool back = false;
            for (const auto& name : components) {
                if (hasComponent(name)) {
                    back = true;
                } else {
                    addExtra(name, back);
                }
                setComponent(name, path[name]);
                exclude.insert(name);
            }
        } else {
            // legacy YAML format does not provide for list of components
            for (const auto& [name, value] : path) {
                if (value.isVector<double>()) {
                    const vector<double>& data = value.asVector<double>();
                    if (data.size() == m_dataSize) {
                        if (!hasComponent(name)) {
                            addExtra(name);
                        }
                        setComponent(name, value);
                        exclude.insert(name);
                    }
                }
            }
        }

        // check that state data are complete
        const auto& nativeState = m_sol->thermo()->nativeState();
        set<string> props;
        set<string> missingProps;
        for (const auto& [name, offset] : nativeState) {
            if (exclude.count(name)) {
                props.insert(name);
            } else {
                missingProps.insert(name);
            }
        }

        set<string> TY = {"T", "Y"};
        if (props == TY && missingProps.count("D") && path.hasKey("pressure")) {
            // missing "D" - Sim1D (Cantera 2.6)
            double P = path["pressure"].asDouble();
            const size_t offset_T = nativeState.find("T")->second;
            const size_t offset_D = nativeState.find("D")->second;
            const size_t offset_Y = nativeState.find("Y")->second;
            for (size_t i = 0; i < m_dataSize; i++) {
                double T = (*m_data)[offset_T + i * m_stride];
                m_sol->thermo()->setState_TPY(
                    T, P, m_data->data() + offset_Y + i * m_stride);
                (*m_data)[offset_D + i * m_stride] = m_sol->thermo()->density();
            }
        } else if (missingProps.size()) {
            throw CanteraError("SolutionArray::readEntry",
                "Incomplete state information: missing '{}'.",
                ba::join(missingProps, "', '"));
        }
    }

    // add meta data
    for (const auto& [name, value] : path) {
        if (!exclude.count(name)) {
            m_meta[name] = value;
        }
    }
}

namespace { // restrict scope of helper functions to local translation unit

template<class T>
AnyValue getSingle(const AnyValue& extra, const vector<int>& slice)
{
    vector<T> data(slice.size());
    const auto& vec = extra.asVector<T>();
    for (size_t k = 0; k < slice.size(); ++k) {
        data[k] = vec[slice[k]];
    }
    AnyValue out;
    out = data;
    return out;
}

template<class T>
AnyValue getMulti(const AnyValue& extra, const vector<int>& slice)
{
    vector<vector<T>> data(slice.size());
    const auto& vec = extra.asVector<vector<T>>();
    for (size_t k = 0; k < slice.size(); ++k) {
        data[k] = vec[slice[k]];
    }
    AnyValue out;
    out = data;
    return out;
}

template<class T>
void setScalar(AnyValue& extra, const AnyValue& data, const vector<int>& slice)
{
    T value = data.as<T>();
    if (extra.isVector<T>()) {
        auto& vec = extra.asVector<T>();
        for (size_t k = 0; k < slice.size(); ++k) {
            vec[slice[k]] = value;
        }
    } else {
        throw CanteraError("SolutionArray::setScalar",
            "Incompatible input data: unable to assign '{}' data to '{}'.",
            data.type_str(), extra.type_str());
    }
}

template<class T>
void setSingle(AnyValue& extra, const AnyValue& data, const vector<int>& slice)
{
    size_t size = slice.size();
    if (extra.vectorSize() == size && data.vectorSize() == size) {
        extra = data; // no slicing necessary; type can change
        return;
    }
    if (extra.matrixShape().first == size && data.vectorSize() == size) {
        extra = data; // no slicing necessary; shape and type can change
        return;
    }
    if (extra.type_str() != data.type_str()) {
        // do not allow changing of data type when slicing
        throw CanteraError("SolutionArray::setSingle",
            "Incompatible types: expected '{}' but received '{}'.",
            extra.type_str(), data.type_str());
    }
    const auto& vData = data.asVector<T>();
    if (vData.size() != size) {
        throw CanteraError("SolutionArray::setSingle",
            "Invalid input data size: expected {} entries but received {}.",
            size, vData.size());
    }
    auto& vec = extra.asVector<T>();
    for (size_t k = 0; k < size; ++k) {
        vec[slice[k]] = vData[k];
    }
}

template<class T>
void setMulti(AnyValue& extra, const AnyValue& data, const vector<int>& slice)
{
    if (!data.isMatrix<T>()) {
        throw CanteraError("SolutionArray::setMulti",
            "Invalid input data shape: inconsistent number of columns.");
    }
    size_t size = slice.size();
    auto [rows, cols] = data.matrixShape();
    if (extra.matrixShape().first == size && rows == size) {
        extra = data; // no slicing necessary; type can change
        return;
    }
    if (extra.vectorSize() == size && rows == size) {
        extra = data; // no slicing necessary; shape and type can change
        return;
    }
    if (extra.type_str() != data.type_str()) {
        // do not allow changing of data type when slicing
        throw CanteraError("SolutionArray::setMulti",
            "Incompatible types: expected '{}' but received '{}'.",
            extra.type_str(), data.type_str());
    }
    if (rows != size) {
        throw CanteraError("SolutionArray::setMulti",
            "Invalid input data shape: expected {} rows but received {}.",
            size, rows);
    }
    if (extra.matrixShape().second != cols) {
        throw CanteraError("SolutionArray::setMulti",
            "Invalid input data shape: expected {} columns but received {}.",
            extra.matrixShape().second, cols);
    }
    const auto& vData = data.asVector<vector<T>>();
    auto& vec = extra.asVector<vector<T>>();
        for (size_t k = 0; k < slice.size(); ++k) {
            vec[slice[k]] = vData[k];
    }
}

template<class T>
void resizeSingle(AnyValue& extra, size_t size, const AnyValue& value)
{
    T defaultValue;
    if (value.is<void>()) {
        defaultValue = vector<T>(1)[0];
    } else {
        defaultValue = value.as<T>();
    }
    extra.asVector<T>().resize(size, defaultValue);
}

template<class T>
void resizeMulti(AnyValue& extra, size_t size, const AnyValue& value)
{
    vector<T> defaultValue;
    if (value.is<void>()) {
        defaultValue = vector<T>(extra.matrixShape().second);
    } else {
        defaultValue = value.as<vector<T>>();
    }
    extra.asVector<vector<T>>().resize(size, defaultValue);
}

template<class T>
void resetSingle(AnyValue& extra, const vector<int>& slice)
{
    T defaultValue = vector<T>(1)[0];
    vector<T>& data = extra.asVector<T>();
    for (size_t k = 0; k < slice.size(); ++k) {
        data[slice[k]] = defaultValue;
    }
}

template<class T>
void resetMulti(AnyValue& extra, const vector<int>& slice)
{
    vector<T> defaultValue = vector<T>(extra.matrixShape().second);
    vector<vector<T>>& data = extra.asVector<vector<T>>();
    for (size_t k = 0; k < slice.size(); ++k) {
        data[slice[k]] = defaultValue;
    }
}

template<class T>
void setAuxiliarySingle(size_t loc, AnyValue& extra, const AnyValue& value)
{
    extra.asVector<T>()[loc] = value.as<T>();
}

template<class T>
void setAuxiliaryMulti(size_t loc, AnyValue& extra, const AnyValue& data)
{
    const auto& value = data.asVector<T>();
    auto& vec = extra.asVector<vector<T>>();
    if (value.size() != vec[loc].size()) {
        throw CanteraError("SolutionArray::setAuxiliaryMulti",
            "New element size {} does not match existing column size {}.",
            value.size(), vec[loc].size());
    }
    vec[loc] = value;
}

} // end unnamed namespace

}
