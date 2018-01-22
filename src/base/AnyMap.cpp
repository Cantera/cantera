//! @file AnyMap.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/AnyMap.h"
#include "cantera/base/yaml.h"
#include "cantera/base/stringUtils.h"

#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;

namespace { // helper functions

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
        if (!node.IsMap()) {
            std::string text = YAML::Dump(node);
            if (text.size() > 300) {
                text.resize(300);
            }
            throw CanteraError("YAML::convert<AnyMap>",
                "YAML node is not a map. Node begins with:\n'''\n{}\n'''", text);
        }
        for (const auto& child : node) {
            std::string key = child.first.as<std::string>();
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
        }
        return false;
    }
};

}

namespace Cantera
{

std::map<std::string, std::string> AnyValue::s_typenames = {
    {typeid(double).name(), "double"},
    {typeid(std::string).name(), "string"},
    {typeid(std::vector<double>).name(), "vector<double>"},
    {typeid(AnyMap).name(), "AnyMap"},
};

// Methods of class AnyValue

AnyValue::AnyValue()
  : m_key()
  , m_value(new boost::any{})
{}

AnyValue::~AnyValue() = default;

AnyValue::AnyValue(AnyValue const& other): m_key(other.m_key),
                                           m_value(new boost::any{*other.m_value}) {
}

AnyValue::AnyValue(AnyValue&& other): m_key(std::move(other.m_key)),
                                      m_value(std::move(other.m_value)) {
}

AnyValue& AnyValue::operator=(AnyValue const& other) {
    if (this == &other)
        return *this;
    m_key = other.m_key;
    m_value.reset(new boost::any{*other.m_value});
    return *this;
}

AnyValue& AnyValue::operator=(AnyValue&& other) {
    if (this == &other)
        return *this;
    m_key = std::move(other.m_key);
    m_value = std::move(other.m_value);
    return *this;
}

AnyValue& AnyValue::operator[](const std::string& key)
{
    return as<AnyMap>()[key];
}

bool AnyValue::hasKey(const std::string& key) const {
    return (is<AnyMap>() && as<AnyMap>().hasKey(key));
}

void AnyValue::setKey(const std::string &key) { m_key = key; }

const std::type_info &AnyValue::type() {
    return m_value->type();
}

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

AnyValue &AnyValue::operator=(double value) {
    *m_value = value;
    return *this;
}

double AnyValue::asDouble() const {
    return as<double>();
}

AnyValue &AnyValue::operator=(bool value) {
    *m_value = value;
    return *this;
}

bool AnyValue::asBool() const {
    return as<bool>();
}

AnyValue &AnyValue::operator=(long int value) {
    *m_value = value;
    return *this;
}

AnyValue &AnyValue::operator=(int value) {
    *m_value = static_cast<long int>(value);
    return *this;
}

long int AnyValue::asInt() const {
    return as<long int>();
}
AnyValue& AnyValue::operator=(const AnyMap& value) {
    *m_value = value;
    return *this;
}

AnyValue& AnyValue::operator=(AnyMap&& value) {
    *m_value = std::move(value);
    return *this;
}

std::string AnyValue::demangle(const std::type_info& type) const
{
    if (s_typenames.find(type.name()) != s_typenames.end()) {
        return s_typenames[type.name()];
    } else {
        return type.name();
    }
}

// Methods of class AnyMap

AnyValue& AnyMap::operator[](const std::string& key)
{
    const auto& slash = boost::ifind_first(key, "/");
    if (!slash) {
        // Simple key
        const auto& iter = m_data.find(key);
        if (iter == m_data.end()) {
            // Create a new key return it
            // NOTE: 'insert' can be replaced with 'emplace' after support for
            // G++ 4.7 is dropped.
            AnyValue& value = m_data.insert({key, AnyValue()}).first->second;
            value.setKey(key);
            return value;
        } else {
            // Return an already-existing item
            return iter->second;
        }
    } else {
        // Split the first slash-delimited part of key and recurse
        std::string head(key.begin(), slash.begin());
        std::string tail(slash.end(), key.end());
        const auto& iter = m_data.find(head);
        if (iter == m_data.end()) {
            // Create a new key
            AnyValue& value = m_data.insert({head, AnyValue()}).first->second;
            value = AnyMap();
            value.setKey(head);
            return value.as<AnyMap>()[tail];
        } else {
            // Return an already existing key
            return iter->second.as<AnyMap>()[tail];
        }
    }
}

const AnyValue& AnyMap::at(const std::string& key) const
{
    const auto& slash = boost::ifind_first(key, "/");
    if (!slash) {
        return m_data.at(key);
    } else {
        std::string head(key.begin(), slash.begin());
        std::string tail(slash.end(), key.end());
        return m_data.at(head).as<AnyMap>().at(tail);
    }
}

bool AnyMap::hasKey(const std::string& key) const
{
    const auto& slash = boost::ifind_first(key, "/");
    if (!slash) {
        return (m_data.find(key) != m_data.end());
    } else {
        std::string head(key.begin(), slash.begin());
        std::string tail(slash.end(), key.end());
        if (m_data.find(head) == m_data.end() || !m_data.at(head).is<AnyMap>()) {
            return false;
        } else {
            return m_data.at(head).as<AnyMap>().hasKey(tail);
        }
    }
}

AnyMap AnyMap::fromYamlString(const std::string& yaml) {
    YAML::Node node = YAML::Load(yaml);
    return node.as<AnyMap>();
}

AnyMap AnyMap::fromYamlFile(const std::string& name) {
    YAML::Node node = YAML::LoadFile(findInputFile(name));
    return node.as<AnyMap>();
}

}
