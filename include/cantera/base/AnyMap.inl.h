//! @file AnyMap.inl.h

#ifndef CT_ANYMAP_INL_H
#define CT_ANYMAP_INL_H

#include "cantera/base/AnyMap.h"

#include <boost/any.hpp>
#include <boost/algorithm/string.hpp>

namespace Cantera
{

// Definitions for AnyValue templated functions

template<class T>
const T &AnyValue::as() const {
    try {
        return boost::any_cast<const T&>(*m_value);
    } catch (boost::bad_any_cast&) {
        if (m_value->type() == typeid(void)) {
            // Values that have not been set are of type 'void'
            throw CanteraError("AnyValue::as", "Key '{}' not found", m_key);
        } else {
            throw CanteraError("AnyValue::as",
                               "Key '{}' contains a '{}',\nnot a '{}'.",
                               m_key, demangle(m_value->type()), demangle(typeid(T)));
        }
    }
}

template<class T>
T &AnyValue::as() {
    try {
        return boost::any_cast<T&>(*m_value);
    } catch (boost::bad_any_cast&) {
        if (m_value->type() == typeid(void)) {
            // Values that have not been set are of type 'void'
            throw CanteraError("AnyValue::as", "Key '{}' not found", m_key);
        } else {
            throw CanteraError("AnyValue::as",
                               "Key '{}' contains a '{}',\nnot a '{}'.",
                               m_key, demangle(m_value->type()), demangle(typeid(T)));
        }
    }
}

template<class T>
bool AnyValue::is() const {
    return m_value->type() == typeid(T);
}

template<class T>
AnyValue &AnyValue::operator=(const std::vector<T> &value) {
    *m_value = value;
    return *this;
}

template<class T>
const std::vector<T> &AnyValue::asVector() const {
    return as<std::vector<T>>();
}

template<class T>
std::vector<T> &AnyValue::asVector() {
    return as<std::vector<T>>();
}

template<class T>
AnyValue& AnyValue::operator=(const std::unordered_map<std::string, T> items) {
    *m_value = AnyMap();
    AnyMap& dest = as<AnyMap>();
    for (const auto& item : items) {
        dest[item.first] = item.second;
    }
    return *this;
}

template<class T>
AnyValue& AnyValue::operator=(const std::map<std::string, T> items) {
    *m_value = AnyMap();
    AnyMap& dest = as<AnyMap>();
    for (const auto& item : items) {
        dest[item.first] = item.second;
    }
    return *this;
}

template<>
inline AnyMap& AnyValue::as<AnyMap>() {
    try {
        // This is where nested AnyMaps are created when the syntax
        // m[key1][key2] is used.
        if (m_value->type() == typeid(void)) {
            *m_value = AnyMap();
        }
        return boost::any_cast<AnyMap&>(*m_value);
    } catch (boost::bad_any_cast&) {
        throw CanteraError("AnyValue::as",
            "value of key '{}' is a '{}',\nnot an 'AnyMap'.",
            m_key, demangle(m_value->type()));
    }
}

template<class T>
std::map<std::string, T> AnyValue::asMap() const
{
    std::map<std::string, T> dest;
    for (const auto& item : as<AnyMap>().m_data) {
        try {
            dest[item.first] = boost::any_cast<T>(*item.second.m_value);
        } catch (boost::bad_any_cast&) {
            throw CanteraError("AnyValue::asMap",
                "Value of key '{}' is not a '{}'",
                item.first, demangle(typeid(T)));
        }
    }
    return dest;
}

}
#endif
