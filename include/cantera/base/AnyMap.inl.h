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
        if (typeid(T) == typeid(double) && m_value->type() == typeid(long int)) {
            // Implicit conversion of long int to double
            *m_value = static_cast<double>(as<long int>());
        }
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
        if (typeid(T) == typeid(double) && m_value->type() == typeid(long int)) {
            // Implicit conversion of long int to double
            *m_value = static_cast<double>(as<long int>());
        }
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
const std::vector<T> &AnyValue::asVector(size_t nMin, size_t nMax) const {
    const auto& v = as<std::vector<T>>();
    checkSize(v, nMin, nMax);
    return v;
}

template<class T>
std::vector<T> &AnyValue::asVector(size_t nMin, size_t nMax) {
    auto& v = as<std::vector<T>>();
    checkSize(v, nMin, nMax);
    return v;
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
        dest[item.first] = item.second.as<T>();
    }
    return dest;
}

template<class T>
void AnyValue::checkSize(const std::vector<T>& v, size_t nMin, size_t nMax) const
{
    if (nMin != npos && nMax == npos && v.size() != nMin) {
        throw CanteraError("AnyValue::checkSize", "Expected array '{}' "
            "to have length {}, but found an array of length {}.",
            m_key, nMin, v.size());
    } else if (nMin != npos && nMax != npos
               && (v.size() < nMin || v.size() > nMax)) {
        throw CanteraError("AnyValue::checkSize",
            "Expected array '{}' to have from {} to {} elements, but found an "
            " array of length {}.", m_key, nMin, nMax, v.size());
    }
}

}
#endif
