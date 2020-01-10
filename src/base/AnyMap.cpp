//! @file AnyMap.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/AnyMap.h"
#include "application.h"
#include "cantera/base/yaml.h"
#include "cantera/base/stringUtils.h"
#ifdef CT_USE_DEMANGLE
  #include <boost/core/demangle.hpp>
#endif

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <mutex>

namespace ba = boost::algorithm;

namespace { // helper functions

std::mutex yaml_cache_mutex;

bool isFloat(const std::string& val)
{
    // This function duplicates the logic of fpValueCheck, but doesn't throw
    // exceptions if the string isn't a float
    std::string str = ba::trim_copy(val);
    if (str.empty()) {
        return false;
    }
    int numDot = 0;
    int numExp = 0;
    int istart = 0;
    char ch = str[0];
    if (ch == '+' || ch == '-') {
        istart = 1;
        if (str.size() == 1) {
            return false;
        }
    }
    for (size_t i = istart; i < str.size(); i++) {
        ch = str[i];
        if (isdigit(ch)) {
        } else if (ch == '.') {
            numDot++;
            if (numDot > 1) {
                return false;
            }
            if (numExp > 0) {
                return false;
            }
        } else if (ch == 'e' || ch == 'E') {
            numExp++;
            if (numExp > 1 || i == str.size() - 1) {
                return false;
            }
            ch = str[i+1];
            if (ch == '+' || ch == '-') {
                if (i + 1 == str.size() - 1) {
                    return false;
                }
                i++;
            }
        } else {
            return false;
        }
    }
    return true;
}

bool isInt(const std::string& val)
{
    std::string str = ba::trim_copy(val);
    if (str.empty()) {
        return false;
    }
    int istart = 0;
    char ch = str[0];
    if (ch == '+' || ch == '-') {
        istart = 1;
        if (str.size() == 1) {
            return false;
        }
    }
    for (size_t i = istart; i < str.size(); i++) {
        if (!isdigit(str[i])) {
            return false;
        }
    }
    return true;
}

bool isBool(const std::string& val) {
    std::string str = ba::trim_copy(val);
    return (val == "true" || val == "True" || val == "false" || val == "False");
}

enum class Type : char {
    Unknown = 0,
    Integer = 1,
    Double = 2,
    String = 4,
    Bool = 8,
    Map = 16,
    Sequence = 32
};

Type operator|(Type lhs, Type rhs)
{
    return Type(static_cast<char>(lhs) | static_cast<char>(rhs));
}

Type elementTypes(const YAML::Node& node)
{
    // See what kinds of elements we have:
    Type types = Type::Unknown;
    for (const auto& el : node) {
        if (el.IsMap()) {
            types = types | Type::Map;
        } else if (el.IsSequence()) {
            types = types | Type::Sequence;
        } else if (el.IsScalar()) {
            std::string nodestr = el.as<std::string>();
            if (isInt(nodestr)) {
                types = types | Type::Integer;
            } else if (isFloat(nodestr)) {
                types = types | Type::Double;
            } else if (isBool(nodestr)) {
                types = types | Type::Bool;
            } else {
                types = types | Type::String;
            }
        }
    }
    return types;
}

} // end anonymous namespace

namespace YAML { // YAML converters

using namespace Cantera;

template<>
struct convert<Cantera::AnyMap> {
    static Node encode(const Cantera::AnyMap& rhs) {
        throw NotImplementedError("AnyMap::encode");
    }
    static bool decode(const Node& node, Cantera::AnyMap& target) {
        target.setLoc(node.Mark().line, node.Mark().column);
        if (node.IsSequence()) {
            // Convert a top-level list to a map with the key "items"
            target["items"] = node.as<AnyValue>();
            return true;
        } else if (!node.IsMap()) {
            std::string text = YAML::Dump(node);
            if (text.size() > 300) {
                text.resize(300);
            }
            throw CanteraError("YAML::convert<AnyMap>",
                "YAML node is not a map. Node begins with:\n'''\n{}\n'''", text);
        }
        for (const auto& child : node) {
            std::string key = child.first.as<std::string>();
            const auto& loc = child.second.Mark();
            target[key].setLoc(loc.line, loc.column);
            if (child.second.IsMap()) {
                target[key] = child.second.as<AnyMap>();
            } else {
                target[key] = child.second.as<AnyValue>();
                target[key].setKey(key);
            }
        }
        return true;
    }
};

template<>
struct convert<Cantera::AnyValue> {
    static Node encode(const Cantera::AnyValue& rhs) {
       throw NotImplementedError("AnyValue::encode");
    }

    static bool decode(const Node& node, Cantera::AnyValue& target) {
        target.setLoc(node.Mark().line, node.Mark().column);
        if (node.IsScalar()) {
            // Scalar nodes are int, doubles, or strings
            std::string nodestr = node.as<std::string>();
            if (isInt(nodestr)) {
                target = intValue(nodestr);
            } else if (isFloat(nodestr)) {
                target = fpValue(nodestr);
            } else if (isBool(nodestr)) {
                target = node.as<bool>();
            } else {
                target = nodestr;
            }
            return true;
        } else if (node.IsSequence()) {
            // Convert sequences of the same element type to vectors of that type
            Type types = elementTypes(node);
            if (types == Type::Integer) {
                target = node.as<std::vector<long int>>();
            } else if (types == (Type::Integer | Type::Double) || types == Type::Double) {
                target = node.as<vector_fp>();
            } else if (types == Type::String) {
                target = node.as<std::vector<std::string>>();
            } else if (types == Type::Bool) {
                target = node.as<std::vector<bool>>();
            } else if (types == Type::Map) {
                target = node.as<std::vector<AnyMap>>();
            } else if (types == Type::Sequence) {
                // Create nested vectors when data types are compatible
                Type subtypes = Type::Unknown;
                for (const auto& el : node) {
                    subtypes = subtypes | elementTypes(el);
                }
                if (subtypes == Type::Integer) {
                    target = node.as<std::vector<std::vector<long int>>>();
                } else if (subtypes == (Type::Integer | Type::Double) || subtypes == Type::Double) {
                    target = node.as<std::vector<std::vector<double>>>();
                } else if (types == Type::String) {
                    target = node.as<std::vector<std::vector<std::string>>>();
                } else if (types == Type::Bool) {
                    target = node.as<std::vector<std::vector<bool>>>();
                } else {
                    target = node.as<std::vector<AnyValue>>();
                }
            } else {
                // If types are different, create a vector of generic values
                target = node.as<std::vector<AnyValue>>();
            }
            return true;
        } else if (node.IsMap()) {
            target = node.as<AnyMap>();
            return true;
        }
        return false;
    }
};

}

namespace Cantera
{

class InputFile
{
public:
    std::string name;
    std::string contents;
};

std::map<std::string, std::string> AnyValue::s_typenames = {
    {typeid(double).name(), "double"},
    {typeid(long int).name(), "long int"},
    {typeid(std::string).name(), "string"},
    {typeid(std::vector<double>).name(), "vector<double>"},
    {typeid(AnyMap).name(), "AnyMap"},
};

std::unordered_map<std::string, std::pair<AnyMap, int>> AnyMap::s_cache;

// Methods of class AnyValue

AnyValue::AnyValue()
  : m_line(-1)
  , m_column(-1)
  , m_key()
  , m_value(new boost::any{})
{}

AnyValue::~AnyValue() = default;

AnyValue::AnyValue(AnyValue const& other)
    : m_line(other.m_line)
    , m_column(other.m_column)
    , m_file(other.m_file)
    , m_key(other.m_key)
    , m_value(new boost::any{*other.m_value})
{
}

AnyValue::AnyValue(AnyValue&& other)
    : m_line(other.m_line)
    , m_column(other.m_column)
    , m_file(std::move(other.m_file))
    , m_key(std::move(other.m_key))
    , m_value(std::move(other.m_value))
{
}

AnyValue& AnyValue::operator=(AnyValue const& other) {
    if (this == &other)
        return *this;
    m_line = other.m_line;
    m_column = other.m_column;
    m_file = other.m_file;
    m_key = other.m_key;
    m_value.reset(new boost::any{*other.m_value});
    return *this;
}

AnyValue& AnyValue::operator=(AnyValue&& other) {
    if (this == &other)
        return *this;
    m_line = other.m_line;
    m_column = other.m_column;
    m_file = std::move(other.m_file);
    m_key = std::move(other.m_key);
    m_value = std::move(other.m_value);
    return *this;
}

AnyValue& AnyValue::operator[](const std::string& key)
{
    return as<AnyMap>()[key];
}

const AnyValue& AnyValue::operator[](const std::string& key) const
{
    return as<AnyMap>().at(key);
}

bool AnyValue::hasKey(const std::string& key) const {
    return (is<AnyMap>() && as<AnyMap>().hasKey(key));
}

void AnyValue::setKey(const std::string &key) { m_key = key; }

const std::type_info &AnyValue::type() const {
    return m_value->type();
}

void AnyValue::setLoc(int line, int column)
{
    m_line = line;
    m_column = column;
}

void AnyValue::setFile(shared_ptr<InputFile>& file)
{
    m_file = file;
    if (is<AnyMap>()) {
        as<AnyMap>().setFile(m_file);
    } else if (is<std::vector<AnyValue>>()) {
        for (auto& item : asVector<AnyValue>()) {
            item.setFile(m_file);
        }
    } else if (is<std::vector<AnyMap>>()) {
        for (auto& item : asVector<AnyMap>()) {
            item.setFile(m_file);
        }
    }
}

std::string AnyValue::type_str() const {
    return demangle(type());
}

bool AnyValue::isScalar() const {
    return is<double>() || is<long int>() || is<std::string>() || is<bool>();
}

AnyValue::AnyValue(const std::string& value) : m_value(new boost::any{value}) {}

AnyValue::AnyValue(const char* value) : m_value(new boost::any{std::string(value)}) {}

AnyValue &AnyValue::operator=(const std::string &value) {
    *m_value = value;
    return *this;
}

AnyValue &AnyValue::operator=(const char *value) {
    *m_value = std::string(value);
    return *this;
}

const std::string &AnyValue::asString() const {
    return as<std::string>();
}

AnyValue::AnyValue(double value) : m_value(new boost::any{value}) {}

AnyValue &AnyValue::operator=(double value) {
    *m_value = value;
    return *this;
}

double& AnyValue::asDouble() {
    return as<double>();
}

const double& AnyValue::asDouble() const {
    return as<double>();
}

AnyValue::AnyValue(bool value) : m_value(new boost::any{value}) {}

AnyValue &AnyValue::operator=(bool value) {
    *m_value = value;
    return *this;
}

bool& AnyValue::asBool() {
    return as<bool>();
}

const bool& AnyValue::asBool() const {
    return as<bool>();
}

AnyValue::AnyValue(long int value) : m_value(new boost::any{value}) {}

AnyValue &AnyValue::operator=(long int value) {
    *m_value = value;
    return *this;
}

AnyValue &AnyValue::operator=(int value) {
    *m_value = static_cast<long int>(value);
    return *this;
}

long int& AnyValue::asInt() {
    return as<long int>();
}

const long int& AnyValue::asInt() const {
    return as<long int>();
}

AnyValue::AnyValue(const AnyMap& value) : m_value(new boost::any{value}) {}

AnyValue& AnyValue::operator=(const AnyMap& value) {
    *m_value = value;
    return *this;
}

AnyValue& AnyValue::operator=(AnyMap&& value) {
    *m_value = std::move(value);
    return *this;
}

std::unordered_map<std::string, const AnyMap*> AnyValue::asMap(
    const std::string& name) const
{
    std::unordered_map<std::string, const AnyMap*> mapped;
    for (const auto& item : asVector<AnyMap>()) {
        auto key = item[name].asString();
        if (mapped.count(key)) {
            throw InputFileError("AnyValue::asMap", *this,
                                 "Duplicate key '{}'", key);
        }
        mapped.emplace(std::make_pair(key, &item));
    }
    return mapped;
}

std::unordered_map<std::string, AnyMap*> AnyValue::asMap(const std::string& name)
{
    std::unordered_map<std::string, AnyMap*> mapped;
    for (auto& item : asVector<AnyMap>()) {
        auto key = item.at(name).asString();
        if (mapped.count(key)) {
            throw InputFileError("AnyValue::asMap", *this,
                                 "Duplicate key '{}'", key);
        }
        mapped.emplace(std::make_pair(key, &item));
    }
    return mapped;
}

AnyMap& AnyValue::getMapWhere(const std::string& key, const std::string& value)
{
    if (value == "") {
        return asVector<AnyMap>().at(0);
    }
    for (auto& item : asVector<AnyMap>()) {
        if (item.hasKey(key) && item[key].asString() == value) {
            return item;
        }
    }
    throw InputFileError("AnyValue::getMapWhere", *this,
        "List does not contain a map where '{}' = '{}'", key, value);
}

void AnyValue::applyUnits(const UnitSystem& units)
{
    if (is<AnyMap>()) {
        // Units declaration applicable to this map
        as<AnyMap>().applyUnits(units);
    } else if (is<std::vector<AnyMap>>()) {
        auto& list = as<std::vector<AnyMap>>();
        if (list.size() && list[0].hasKey("units") && list[0].size() == 1) {
            // First item in the list is a units declaration, which applies to
            // the items in the list
            UnitSystem newUnits = units;
            newUnits.setDefaults(list[0]["units"].asMap<std::string>());
            list[0].m_data.erase("units");
            for (auto& item : list) {
                // Any additional units declarations are errors
                if (item.size() == 1 && item.hasKey("units")) {
                    throw InputFileError("AnyValue::applyUnits", item,
                        "Found units entry as not the first item in a list.");
                }
                item.applyUnits(newUnits);
            }
            // Remove the "units" map after it has been applied
            list.erase(list.begin());
        } else {
            // Simple downward propagation of the current units
            for (auto& item : list) {
                // Any later units declarations are errors
                if (item.size() == 1 && item.hasKey("units")) {
                    throw InputFileError("AnyValue::applyUnits", item,
                        "Found units entry as not the first item in a list.");
                }
                item.applyUnits(units);
            }
        }
    }

}

std::string AnyValue::demangle(const std::type_info& type) const
{
    if (s_typenames.find(type.name()) != s_typenames.end()) {
        return s_typenames[type.name()];
    } else {
        #ifdef CT_USE_DEMANGLE
            return boost::core::demangle(type.name());
        #else
            return type.name();
        #endif
    }
}

// Explicit template specializations to allow certain conversions

template<>
const std::vector<AnyValue>& AnyValue::asVector<AnyValue>(size_t nMin, size_t nMax) const
{
    if (!is<std::vector<AnyValue>>()) {
        std::vector<AnyValue> v;
        if (is<std::vector<double>>()) {
            for (const auto& el : asVector<double>()) {
                v.push_back(AnyValue(el));
            }
            *m_value = v;
        } else if (is<std::vector<long int>>()) {
            for (const auto& el : asVector<long int>()) {
                v.push_back(AnyValue(el));
            }
            *m_value = v;
        } else if (is<std::vector<std::string>>()) {
            for (const auto& el : asVector<std::string>()) {
                v.push_back(AnyValue(el));
            }
            *m_value = v;
        }
        // If none of these special cases match, the value won't be replaced,
        // and an exception will be thrown.
    }
    const auto& vv = as<std::vector<AnyValue>>();
    checkSize(vv, nMin, nMax);
    return vv;
}

template<>
std::vector<AnyValue>& AnyValue::asVector<AnyValue>(size_t nMin, size_t nMax)
{
    auto& v = const_cast<std::vector<AnyValue>&>(
        const_cast<const AnyValue*>(this)->asVector<AnyValue>());
    checkSize(v, nMin, nMax);
    return v;
}

template<>
const std::vector<double>& AnyValue::asVector<double>(size_t nMin, size_t nMax) const
{
    if (is<std::vector<long int>>()) {
        std::vector<double> v;
        for (const auto& el : asVector<long int>()) {
            v.push_back(el);
        }
        *m_value = v;
    }
    const auto& vv = as<std::vector<double>>();
    checkSize(vv, nMin, nMax);
    return vv;
}

template<>
std::vector<double>& AnyValue::asVector<double>(size_t nMin, size_t nMax)
{
    if (is<std::vector<long int>>()) {
        std::vector<double> v;
        for (const auto& el : asVector<long int>()) {
            v.push_back(el);
        }
        *m_value = v;
    }
    auto& vv = as<std::vector<double>>();
    checkSize(vv, nMin, nMax);
    return vv;
}

template<>
const std::vector<vector_fp>& AnyValue::asVector<vector_fp>(size_t nMin, size_t nMax) const
{
    if (is<std::vector<std::vector<long int>>>()) {
        std::vector<vector_fp> v;
        for (const auto& outer : asVector<std::vector<long int>>()) {
            v.push_back(vector_fp());
            for (const auto& inner : outer) {
                v.back().push_back(inner);
            }
        }
        *m_value = v;
    }
    const auto& vv = as<std::vector<vector_fp>>();
    checkSize(vv, nMin, nMax);
    return vv;
}

template<>
std::vector<vector_fp>& AnyValue::asVector<vector_fp>(size_t nMin, size_t nMax)
{
    if (is<std::vector<std::vector<long int>>>()) {
        std::vector<vector_fp> v;
        for (const auto& outer : asVector<std::vector<long int>>()) {
            v.push_back(vector_fp());
            for (const auto& inner : outer) {
                v.back().push_back(inner);
            }
        }
        *m_value = v;
    }
    auto& vv = as<std::vector<vector_fp>>();
    checkSize(vv, nMin, nMax);
    return vv;
}


// Methods of class AnyMap

AnyValue& AnyMap::operator[](const std::string& key)
{
    const auto& iter = m_data.find(key);
    if (iter == m_data.end()) {
        // Create a new key return it
        // NOTE: 'insert' can be replaced with 'emplace' after support for
        // G++ 4.7 is dropped.
        AnyValue& value = m_data.insert({key, AnyValue()}).first->second;
        value.setKey(key);
        if (m_file) {
            // Approximate location, useful mainly if this insertion is going to
            // immediately result in an error that needs to be reported.
            value.setLoc(m_line, m_column);
            value.setFile(m_file);
        }
        return value;
    } else {
        // Return an already-existing item
        return iter->second;
    }
}

const AnyValue& AnyMap::operator[](const std::string& key) const
{
    try {
        return m_data.at(key);
    } catch (std::out_of_range&) {
        throw InputFileError("AnyMap::operator[]", *this,
            "Key '{}' not found.\nExisting keys: {}", key, keys_str());
    }
}

const AnyValue& AnyMap::at(const std::string& key) const
{
    try {
        return m_data.at(key);
    } catch (std::out_of_range&) {
        throw InputFileError("AnyMap::at", *this,
            "Key '{}' not found.\nExisting keys: {}", key, keys_str());
    }
}

bool AnyMap::hasKey(const std::string& key) const
{
    return (m_data.find(key) != m_data.end());
}

void AnyMap::erase(const std::string& key)
{
    m_data.erase(key);
}

std::string AnyMap::keys_str() const
{
    fmt::memory_buffer b;
    auto iter = this->begin();
    if (iter != this->end()) {
        format_to(b, "{}", iter->first);
        ++iter;
    }
    while (iter != this->end()) {
        format_to(b, ", {}", iter->first);
        ++iter;
    }
    return to_string(b);
}

void AnyMap::setLoc(int line, int column)
{
    m_line = line;
    m_column = column;
}

void AnyMap::setFile(shared_ptr<InputFile>& file)
{
    m_file = file;
    for (auto& item : m_data) {
        item.second.setFile(m_file);
    }
}

void AnyMap::setFileName(const std::string& filename)
{
    auto info = make_shared<InputFile>();
    info->name = filename;
    setFile(info);
}

void AnyMap::setFileContents(const std::string& contents)
{
    auto info = make_shared<InputFile>();
    info->contents = contents;
    setFile(info);
}

bool AnyMap::getBool(const std::string& key, bool default_) const
{
    return (hasKey(key)) ? m_data.at(key).asBool() : default_;
}

double AnyMap::getDouble(const std::string& key, double default_) const
{
    return (hasKey(key)) ? m_data.at(key).asDouble() : default_;
}

long int AnyMap::getInt(const std::string& key, long int default_) const
{
    return (hasKey(key)) ? m_data.at(key).asInt() : default_;
}

const std::string& AnyMap::getString(const std::string& key,
                                     const std::string& default_) const
{
    return (hasKey(key)) ? m_data.at(key).asString() : default_;
}

double AnyMap::convert(const std::string& key, const std::string& dest) const
{
    return units().convert(at(key), dest);
}

double AnyMap::convert(const std::string& key, const Units& dest) const
{
    return units().convert(at(key), dest);
}

double AnyMap::convert(const std::string& key, const std::string& dest,
                       double default_) const
{
    if (hasKey(key)) {
        return units().convert(at(key), dest);
    } else {
        return default_;
    }
}

vector_fp AnyMap::convertVector(const std::string& key, const std::string& dest,
                                size_t nMin, size_t nMax) const
{
    return units().convert(at(key).asVector<AnyValue>(nMin, nMax), dest);
}

void AnyMap::applyUnits(const UnitSystem& units) {
    m_units = units;

    if (hasKey("units")) {
        m_units.setDefaults(at("units").asMap<std::string>());
        m_data.erase("units");
    }
    for (auto& item : m_data) {
        item.second.applyUnits(m_units);
    }
}

AnyMap AnyMap::fromYamlString(const std::string& yaml) {
    AnyMap amap;
    try {
        YAML::Node node = YAML::Load(yaml);
        amap = node.as<AnyMap>();
    } catch (YAML::Exception& err) {
        AnyMap fake;
        fake.setLoc(err.mark.line, err.mark.column);
        fake.setFileContents(yaml);
        throw InputFileError("AnyMap::fromYamlString", fake, err.msg);
    }
    amap.setFileContents(yaml);
    amap.applyUnits(UnitSystem());
    return amap;
}

AnyMap AnyMap::fromYamlFile(const std::string& name,
                            const std::string& parent_name)
{
    std::string fullName;
    // See if a file with this name exists in a path relative to the parent file
    size_t islash = parent_name.find_last_of("/\\");
    if (islash != npos) {
        std::string parent_path = parent_name.substr(0, islash);
        if (std::ifstream(parent_path + "/" + name).good()) {
            fullName = parent_path + "/" + name;
        }
    }
    // Otherwise, search the Cantera include path for the file
    if (fullName.empty()) {
        fullName = findInputFile(name);
    }

    // Check for an already-parsed YAML file with the same last-modified time,
    // and return that if possible
    int mtime = get_modified_time(fullName);
    std::unique_lock<std::mutex> lock(yaml_cache_mutex);
    auto iter = s_cache.find(fullName);
    if (iter != s_cache.end() && iter->second.second == mtime) {
        return iter->second.first;
    }

    if (!std::ifstream(fullName).good()) {
        throw CanteraError("AnyMap::fromYamlFile", "Input file '{}' not found "
            "on the Cantera search path.", name);
    }

    // Generate an AnyMap from the YAML file and store it in the cache
    auto& cache_item = s_cache[fullName];
    cache_item.second = mtime;
    try {
        YAML::Node node = YAML::LoadFile(fullName);
        cache_item.first = node.as<AnyMap>();
        cache_item.first.setFileName(fullName);
        cache_item.first.applyUnits(UnitSystem());
    } catch (YAML::Exception& err) {
        s_cache.erase(fullName);
        AnyMap fake;
        fake.setLoc(err.mark.line, err.mark.column);
        fake.setFileName(fullName);
        throw InputFileError("AnyMap::fromYamlFile", fake, err.msg);
    } catch (CanteraError&) {
        s_cache.erase(fullName);
        throw;
    }
    cache_item.first["__file__"] = fullName;

    if (cache_item.first.hasKey("deprecated")) {
        warn_deprecated(fullName, cache_item.first["deprecated"].asString());
    }

    // Return a copy of the AnyMap
    return cache_item.first;
}

AnyMap::const_iterator begin(const AnyValue& v) {
    return v.as<AnyMap>().begin();
}

AnyMap::const_iterator end(const AnyValue& v) {
    return v.as<AnyMap>().end();
}

std::string InputFileError::formatError(const std::string& message,
                                        int lineno, int column,
                                        const shared_ptr<InputFile>& file)
{
    if (!file) {
        return message;
    }

    fmt::memory_buffer b;
    format_to(b, "Error on line {} of", lineno+1);
    if (file->name.empty()) {
        format_to(b, " input string:\n");
    } else {
        format_to(b, " {}:\n", file->name);
    }
    format_to(b, "{}\n", message);
    format_to(b, "|  Line |\n");
    if (file->contents.empty()) {
        std::ifstream infile(findInputFile(file->name));
        std::stringstream buffer;
        buffer << infile.rdbuf();
        file->contents = buffer.str();
    }
    std::string line;
    int i = 0;
    std::stringstream contents(file->contents);
    while (std::getline(contents, line)) {
        if (lineno == i) {
            format_to(b, "> {: 5d} > {}\n", i+1, line);
            format_to(b, "{:>{}}\n", "^", column + 11);
        } else if (lineno + 4 > i && lineno < i + 6) {
            format_to(b, "| {: 5d} | {}\n", i+1, line);
        }
        i++;
    }
    return to_string(b);
}

}
