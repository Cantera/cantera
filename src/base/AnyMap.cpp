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
#include <unordered_set>

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

Cantera::AnyValue Empty;

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
                try {
                    target = node.as<long int>();
                } catch (YAML::BadConversion&) {
                    // This exception is raised if the value doesn't fit in a
                    // long int, in which case we would rather store it
                    // (possibly inexactly) as a double.
                    target = node.as<double>();
                }
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

std::map<std::string, std::string> AnyValue::s_typenames = {
    {typeid(double).name(), "double"},
    {typeid(long int).name(), "long int"},
    {typeid(std::string).name(), "string"},
    {typeid(std::vector<double>).name(), "vector<double>"},
    {typeid(AnyMap).name(), "AnyMap"},
};

std::unordered_map<std::string, std::pair<AnyMap, int>> AnyMap::s_cache;

// Methods of class AnyBase

AnyBase::AnyBase()
    : m_line(-1)
    , m_column(-1)
{}

void AnyBase::setLoc(int line, int column)
{
    m_line = line;
    m_column = column;
}

const AnyValue& AnyBase::getMetadata(const std::string& key) const
{
    if (m_metadata && m_metadata->hasKey(key)) {
        return m_metadata->at(key);
    } else {
        return Empty;
    }
}

// Methods of class AnyValue

AnyValue::AnyValue()
  : m_key()
  , m_value(new boost::any{})
  , m_equals(eq_comparer<size_t>)
{}

AnyValue::~AnyValue() = default;

AnyValue::AnyValue(AnyValue const& other)
    : AnyBase(other)
    , m_key(other.m_key)
    , m_value(new boost::any{*other.m_value})
    , m_equals(other.m_equals)
{
}

AnyValue::AnyValue(AnyValue&& other)
    : AnyBase(std::move(other))
    , m_key(std::move(other.m_key))
    , m_value(std::move(other.m_value))
    , m_equals(std::move(other.m_equals))
{
}

AnyValue& AnyValue::operator=(AnyValue const& other) {
    if (this == &other) {
        return *this;
    }
    AnyBase::operator=(*this);
    m_key = other.m_key;
    m_value.reset(new boost::any{*other.m_value});
    m_equals = other.m_equals;
    return *this;
}

AnyValue& AnyValue::operator=(AnyValue&& other) {
    if (this == &other) {
        return *this;
    }
    AnyBase::operator=(std::move(other));
    m_key = std::move(other.m_key);
    m_value = std::move(other.m_value);
    m_equals = std::move(other.m_equals);
    return *this;
}

bool AnyValue::operator==(const AnyValue& other) const
{
    return m_equals(*m_value, *other.m_value);
}

bool AnyValue::operator!=(const AnyValue& other) const
{
    return !m_equals(*m_value, *other.m_value);
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

void AnyValue::propagateMetadata(shared_ptr<AnyMap>& metadata)
{
    m_metadata = metadata;
    if (is<AnyMap>()) {
        as<AnyMap>().propagateMetadata(m_metadata);
    } else if (is<std::vector<AnyValue>>()) {
        for (auto& item : asVector<AnyValue>()) {
            item.propagateMetadata(m_metadata);
        }
    } else if (is<std::vector<AnyMap>>()) {
        for (auto& item : asVector<AnyMap>()) {
            item.propagateMetadata(m_metadata);
        }
    }
}

std::string AnyValue::type_str() const {
    return demangle(type());
}

bool AnyValue::isScalar() const {
    return is<double>() || is<long int>() || is<std::string>() || is<bool>();
}

// Specializations for "std::string" and "const char*"

AnyValue::AnyValue(const std::string& value)
    : m_value(new boost::any{value})
    , m_equals(eq_comparer<std::string>)
{}

AnyValue::AnyValue(const char* value)
    : m_value(new boost::any{std::string(value)})
    , m_equals(eq_comparer<std::string>)
{}

AnyValue &AnyValue::operator=(const std::string &value) {
    *m_value = value;
    m_equals = eq_comparer<std::string>;
    return *this;
}

AnyValue &AnyValue::operator=(const char *value) {
    *m_value = std::string(value);
    m_equals = eq_comparer<std::string>;
    return *this;
}

const std::string &AnyValue::asString() const {
    return as<std::string>();
}

bool AnyValue::operator==(const std::string& other) const
{
    if (m_value->type() == typeid(std::string)) {
        return boost::any_cast<std::string>(*m_value) == other;
    } else {
        return false;
    }
}

bool AnyValue::operator!=(const std::string& other) const
{
    return !(*this == other);
}

bool operator==(const std::string& lhs, const AnyValue& rhs)
{
    return rhs == lhs;
}

bool operator!=(const std::string& lhs, const AnyValue& rhs)
{
    return rhs != lhs;
}

// Specializations for "double"

AnyValue::AnyValue(double value)
    : m_value(new boost::any{value})
    , m_equals(eq_comparer<double>)
{}

AnyValue &AnyValue::operator=(double value) {
    *m_value = value;
    m_equals = eq_comparer<double>;
    return *this;
}

double& AnyValue::asDouble() {
    return as<double>();
}

const double& AnyValue::asDouble() const {
    return as<double>();
}

bool AnyValue::operator==(const double& other) const
{
    if (m_value->type() == typeid(double)) {
        return boost::any_cast<double>(*m_value) == other;
    } else if (m_value->type() == typeid(long int)) {
        return boost::any_cast<long int>(*m_value) == other;
    } else {
        return false;
    }
}

bool AnyValue::operator!=(const double& other) const
{
    return !(*this == other);
}

bool operator==(const double& lhs, const AnyValue& rhs)
{
    return rhs == lhs;
}

bool operator!=(const double& lhs, const AnyValue& rhs)
{
    return rhs != lhs;
}

// Specializations for "bool"

AnyValue::AnyValue(bool value)
    : m_value(new boost::any{value})
    , m_equals(eq_comparer<bool>)
{}

AnyValue &AnyValue::operator=(bool value) {
    *m_value = value;
    m_equals = eq_comparer<bool>;
    return *this;
}

bool& AnyValue::asBool() {
    return as<bool>();
}

const bool& AnyValue::asBool() const {
    return as<bool>();
}

// Specializations for "long int" and "int"

AnyValue::AnyValue(long int value)
    : m_value(new boost::any{value})
    , m_equals(eq_comparer<long int>)
{}

AnyValue::AnyValue(int value)
    : m_value(new boost::any{static_cast<long int>(value)})
    , m_equals(eq_comparer<long int>)
{}

AnyValue &AnyValue::operator=(long int value) {
    *m_value = value;
    m_equals = eq_comparer<long int>;
    return *this;
}

AnyValue &AnyValue::operator=(int value) {
    *m_value = static_cast<long int>(value);
    m_equals = eq_comparer<long int>;
    return *this;
}

long int& AnyValue::asInt() {
    return as<long int>();
}

const long int& AnyValue::asInt() const {
    return as<long int>();
}

bool AnyValue::operator==(const long int& other) const
{
    if (m_value->type() == typeid(long int)) {
        return boost::any_cast<long int>(*m_value) == other;
    } else if (m_value->type() == typeid(double)) {
        return boost::any_cast<double>(*m_value) == other;
    } else {
        return false;
    }
}

bool AnyValue::operator!=(const long int& other) const
{
    return !(*this == other);
}

bool AnyValue::operator==(const int& other) const
{
    return *this == static_cast<long int>(other);
}

bool AnyValue::operator!=(const int& other) const
{
    return *this != static_cast<long int>(other);
}

bool operator==(const long int& lhs, const AnyValue& rhs)
{
    return rhs == lhs;
}

bool operator!=(const long int& lhs, const AnyValue& rhs)
{
    return rhs != lhs;
}

bool operator==(const int& lhs, const AnyValue& rhs)
{
    return rhs == lhs;
}

bool operator!=(const int& lhs, const AnyValue& rhs)
{
    return rhs != lhs;
}

// Specializations for "AnyMap"

AnyValue::AnyValue(const AnyMap& value)
    : m_value(new boost::any{value})
    , m_equals(eq_comparer<AnyMap>)
{}

AnyValue& AnyValue::operator=(const AnyMap& value) {
    *m_value = value;
    m_equals = eq_comparer<AnyMap>;
    return *this;
}

AnyValue& AnyValue::operator=(AnyMap&& value) {
    *m_value = std::move(value);
    m_equals = eq_comparer<AnyMap>;
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

const AnyMap& AnyValue::getMapWhere(const std::string& key, const std::string& value) const
{
    if (is<std::vector<AnyMap>>()) {
        if (value == "") {
            return asVector<AnyMap>().at(0);
        }
        for (auto& item : asVector<AnyMap>()) {
            if (item.hasKey(key) && item[key] == value) {
                return item;
            }
        }
        throw InputFileError("AnyValue::getMapWhere", *this,
            "List does not contain a map where '{}' = '{}'", key, value);
    } else if (is<AnyMap>()) {
        if (value == "" || (hasKey(key) && as<AnyMap>()[key] == value)) {
            return as<AnyMap>();
        } else {
            throw InputFileError("AnyValue::getMapWhere", *this,
                "Map does not contain a key where '{}' = '{}'", key, value);
        }
    } else if (is<void>()) {
        throw InputFileError("AnyValue::getMapWhere", *this,
            "Key '{}' not found", m_key);
    } else {
        throw InputFileError("AnyValue::getMapWhere", *this,
            "Element is not a mapping or list of mappings");
    }
}

AnyMap& AnyValue::getMapWhere(const std::string& key, const std::string& value,
                              bool create)
{
    if (is<std::vector<AnyMap>>()) {
        if (value == "") {
            return asVector<AnyMap>().at(0);
        }
        for (auto& item : asVector<AnyMap>()) {
            if (item.hasKey(key) && item[key] == value) {
                return item;
            }
        }
        if (create) {
            // If the map wasn't found, insert it
            auto& vec = asVector<AnyMap>();
            AnyMap child;
            child[key] = value;
            vec.push_back(std::move(child));
            return vec.back();
        } else {
            throw InputFileError("AnyValue::getMapWhere", *this,
                "List does not contain a map where '{}' = '{}'", key, value);
        }
    } else if (is<AnyMap>()) {
        if (value == "" || (hasKey(key) && as<AnyMap>()[key] == value)) {
            return as<AnyMap>();
        } else if (create) {
            AnyMap newChild;
            newChild[key] = value;
            std::vector<AnyMap> nodes{std::move(as<AnyMap>()), std::move(newChild)};
            operator=(std::move(nodes));
            return asVector<AnyMap>().back();
        } else {
            throw InputFileError("AnyValue::getMapWhere", *this,
                "Map does not contain a key where '{}' = '{}'", key, value);
        }
    } else if (is<void>() && create) {
        AnyMap child;
        child[key] = value;
        operator=(std::move(child));
        return as<AnyMap>();
    } else if (is<void>()) {
        throw InputFileError("AnyValue::getMapWhere", *this,
            "Key '{}' not found", m_key);
    } else {
        throw InputFileError("AnyValue::getMapWhere", *this,
            "Element is not a mapping or list of mappings");
    }
}

bool AnyValue::hasMapWhere(const std::string& key, const std::string& value) const
{
    if (is<std::vector<AnyMap>>()) {
        if (value == "") {
            return true;
        }
        for (auto& item : asVector<AnyMap>()) {
            if (item.hasKey(key) && item[key] == value) {
                return true;
            }
        }
        return false;
    } else if (is<AnyMap>()) {
        if (value == "" || (hasKey(key) && as<AnyMap>()[key] == value)) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
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
    m_equals = eq_comparer<std::vector<AnyValue>>;
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
    m_equals = eq_comparer<std::vector<double>>;
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
    m_equals = eq_comparer<std::vector<double>>;
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
    m_equals = eq_comparer<std::vector<vector_fp>>;
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
    m_equals = eq_comparer<std::vector<vector_fp>>;
    checkSize(vv, nMin, nMax);
    return vv;
}

template<>
const std::vector<AnyMap>& AnyValue::asVector<AnyMap>(size_t nMin, size_t nMax) const
{
    if (is<AnyMap>()) {
        std::vector<AnyMap> v;
        v.push_back(std::move(as<AnyMap>()));
        *m_value = std::move(v);
    } else if (is<std::vector<AnyValue>>() && asVector<AnyValue>().empty()) {
        *m_value = std::vector<AnyMap>();
    }
    const auto& vv = as<std::vector<AnyMap>>();
    checkSize(vv, nMin, nMax);
    return vv;
}

template<>
std::vector<AnyMap>& AnyValue::asVector<AnyMap>(size_t nMin, size_t nMax)
{
    if (is<AnyMap>()) {
        std::vector<AnyMap> v;
        v.push_back(std::move(as<AnyMap>()));
        *m_value = std::move(v);
    } else if (is<std::vector<AnyValue>>() && asVector<AnyValue>().empty()) {
        *m_value = std::vector<AnyMap>();
    }
    auto& vv = as<std::vector<AnyMap>>();
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
        if (m_metadata) {
            // Approximate location, useful mainly if this insertion is going to
            // immediately result in an error that needs to be reported.
            value.setLoc(m_line, m_column);
            value.propagateMetadata(m_metadata);
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

void AnyMap::clear()
{
    m_data.clear();
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

void AnyMap::propagateMetadata(shared_ptr<AnyMap>& metadata)
{
    m_metadata = metadata;
    for (auto& item : m_data) {
        item.second.propagateMetadata(m_metadata);
    }
}

void AnyMap::setMetadata(const std::string& key, const AnyValue& value)
{
    if (m_metadata) {
        // Fork the metadata tree at this point to avoid affecting parent nodes
        m_metadata = make_shared<AnyMap>(*m_metadata);
    } else {
        m_metadata = make_shared<AnyMap>();
    }
    (*m_metadata)[key] = value;
    propagateMetadata(m_metadata);
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

AnyMap::Iterator::Iterator(
    const std::unordered_map<std::string, AnyValue>::const_iterator& start,
    const std::unordered_map<std::string, AnyValue>::const_iterator& stop)
{
    m_iter = start;
    m_stop = stop;
    while (m_iter != m_stop
           && ba::starts_with(m_iter->first, "__")
           && ba::ends_with(m_iter->first, "__")) {
        ++m_iter;
    }
}

AnyMap::Iterator& AnyMap::Iterator::operator++()
{
    ++m_iter;
    while (m_iter != m_stop
           && ba::starts_with(m_iter->first, "__")
           && ba::ends_with(m_iter->first, "__")) {
        ++m_iter;
    }
    return *this;
}

bool AnyMap::operator==(const AnyMap& other) const
{
    // First, make sure that 'other' has all of the non-hidden keys that are in
    // this map
    for (auto& item : *this) {
        if (!other.hasKey(item.first)) {
            return false;
        }
    }
    // Then check for equality, using the non-hidden keys from 'other'
    for (auto & item : other) {
        if (!hasKey(item.first) || item.second != at(item.first)) {
            return false;
        }
    }
    return true;
}

bool AnyMap::operator!=(const AnyMap& other) const
{
    return m_data != other.m_data;
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
        fake.setMetadata("file-contents", AnyValue(yaml));
        throw InputFileError("AnyMap::fromYamlString", fake, err.msg);
    }
    amap.setMetadata("file-contents", AnyValue(yaml));
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
        cache_item.first.setMetadata("filename", AnyValue(fullName));
        cache_item.first.applyUnits(UnitSystem());
    } catch (YAML::Exception& err) {
        s_cache.erase(fullName);
        AnyMap fake;
        fake.setLoc(err.mark.line, err.mark.column);
        fake.setMetadata("filename", AnyValue(fullName));
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

AnyMap::Iterator begin(const AnyValue& v) {
    return v.as<AnyMap>().begin();
}

AnyMap::Iterator end(const AnyValue& v) {
    return v.as<AnyMap>().end();
}

namespace {
void formatInputFile(fmt::memory_buffer& b, const shared_ptr<AnyMap>& metadata,
        const std::string& filename, int lineno, int column, int lineno2=-1, int column2=-1)
{
    if (lineno2 == -1) {
        lineno2 = lineno;
        column2 = column;
    }

    format_to(b, "|  Line |\n");
    if (!metadata->hasKey("file-contents")) {
        std::ifstream infile(findInputFile(filename));
        std::stringstream buffer;
        buffer << infile.rdbuf();
        (*metadata)["file-contents"] = buffer.str();
    }
    std::string line;
    int i = 0;
    int lastShown = -1;
    std::stringstream contents((*metadata)["file-contents"].asString());
    while (std::getline(contents, line)) {
        if (i == lineno || i == lineno2) {
            format_to(b, "> {: 5d} > {}\n", i+1, line);
            format_to(b, "{:>{}}\n", "^", column + 11);
            lastShown = i;
        } else if ((lineno + 4 > i && lineno < i + 6) ||
                   (lineno2 + 4 > i && lineno2 < i + 6)) {
            if (lastShown >= 0 && i - lastShown > 1) {
                format_to(b, "...\n");
            }
            format_to(b, "| {: 5d} | {}\n", i+1, line);
            lastShown = i;
        }
        i++;
    }
}
}

std::string InputFileError::formatError(const std::string& message,
                                        int lineno, int column,
                                        const shared_ptr<AnyMap>& metadata)
{
    if (!metadata) {
        return message;
    }
    std::string filename = metadata->getString("filename", "input string");

    fmt::memory_buffer b;
    format_to(b, "Error on line {} of {}:\n{}\n", lineno+1, filename, message);
    formatInputFile(b, metadata, filename, lineno, column);
    return to_string(b);
}

std::string InputFileError::formatError2(const std::string& message,
                                         int line1, int column1,
                                         const shared_ptr<AnyMap>& metadata1,
                                         int line2, int column2,
                                         const shared_ptr<AnyMap>& metadata2)
{
    if (!metadata1 || !metadata2) {
        return message;
    }
    std::string filename1 = metadata1->getString("filename", "input string");
    std::string filename2 = metadata2->getString("filename", "input string");

    fmt::memory_buffer b;
    if (filename1 == filename2) {
        format_to(b, "Error on lines {} and {} of {}:\n",
                  std::min(line1, line2) + 1, std::max(line1, line2) + 1,
                  filename1);
        format_to(b, "{}\n", message);
        formatInputFile(b, metadata1, filename1, line1, column1, line2, column2);
    } else {
        format_to(b, "Error on line {} of {} and line {} of {}:\n{}\n",
                  line1+1, filename1, line2+1, filename2, message);
        formatInputFile(b, metadata1, filename1, line1, column1);
        format_to(b, "\n");
        formatInputFile(b, metadata2, filename2, line2, column2);
    }

    return to_string(b);
}

}
