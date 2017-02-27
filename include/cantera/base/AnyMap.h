//! @file AnyMap.h

#ifndef CT_ANYMAP_H
#define CT_ANYMAP_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

#include <boost/any.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <unordered_map>

namespace Cantera
{

class AnyMap;

//! A wrapper for a variable whose type is determined at runtime
/*!
 * Instances of AnyValue are used as values in an AnyMap. Values are converted
 * to a concrete type using the templated as() method or convenience methods
 * such as asString() and asDouble(). See AnyMap for usage examples.
 *
 * Elements are set using assignment, and the assignment operator has been
 * overloaded for specific types so that only those types are allowed to be
 * used in an AnyValue.
 */
class AnyValue
{
public:
    AnyValue& operator[](const std::string& key);

    bool hasKey(const std::string& key) const;

    // The value knows the name of its corresponding key in order to provide
    // comprehensible error messages.
    void setKey(const std::string& key) { m_key = key; };

    template<class T>
    const T& as() const {
        try {
            return boost::any_cast<const T&>(m_value);
        } catch (boost::bad_any_cast&) {
            if (m_value.type() == typeid(void)) {
                // Values that have not been set are of type 'void'
                throw CanteraError("AnyValue::as", "Key '{}' not found", m_key);
            } else {
                throw CanteraError("AnyValue::as",
                    "Key '{}' contains a '{}',\nnot a '{}'.",
                    m_key, demangle(m_value.type()), demangle(typeid(T)));
            }
        }
    }

    template<class T>
    T& as() {
        try {
            return boost::any_cast<T&>(m_value);
        } catch (boost::bad_any_cast&) {
            if (m_value.type() == typeid(void)) {
                // Values that have not been set are of type 'void'
                throw CanteraError("AnyValue::as", "Key '{}' not found", m_key);
            } else {
                throw CanteraError("AnyValue::as",
                    "Key '{}' contains a '{}',\nnot a '{}'.",
                    m_key, demangle(m_value.type()), demangle(typeid(T)));
            }
        }
    }

    const std::type_info& type() {
        return m_value.type();
    }

    template<class T>
    bool is() const {
        return m_value.type() == typeid(T);
    }

    AnyValue& operator=(const std::string& value) {
        m_value = value;
        return *this;
    }
    const std::string& asString() const {
        return as<std::string>();
    }

    AnyValue& operator=(double value) {
        m_value = value;
        return *this;
    }
    double asDouble() const {
        return as<double>();
    }

    template<class T>
    AnyValue& operator=(const std::vector<T>& value) {
        m_value = value;
        return *this;
    }
    template<class T>
    const std::vector<T>& asVector() const {
        return as<std::vector<T>>();
    }
    template<class T>
    std::vector<T>& asVector() {
        return as<std::vector<T>>();
    }

    AnyValue& operator=(const AnyMap& value);
    AnyValue& operator=(AnyMap&& value);

    template<class T>
    AnyValue& operator=(const std::unordered_map<std::string, T> items);

    template<class T>
    AnyValue& operator=(const std::map<std::string, T> items);

    template<class T>
    std::map<std::string, T> asMap();

private:
    std::string demangle(const std::type_info& type) const;

    std::string m_key;
    boost::any m_value;
    static std::map<std::string, std::string> s_typenames;
};

//! A map of string keys to values whose type can vary at runtime
/*!
 * Values in an AnyMap are held by instances of AnyValue. Instances of AnyMap
 * can be nested to form a tree.
 *
 * ## Setting elements
 *
 * ```
 * AnyMap breakfast;
 * breakfast["spam"] = 123.4; // Creates a value of type 'double'
 * breakfast["eggs"] = "scrambled"; // Creates a value of type 'std::string'
 *
 * // Create a nested AnyMap named "beans" which has a key named "baked"
 * // whose value is a vector<double>
 * std::vector<double> v{3.14, 1.59, 2.65};
 * breakfast["beans/baked"] = v;
 * // Equivalently:
 * breakfast["beans"]["baked"] = v;
 *
 * // Create a nested AnyMap with values of the same type
 * std::map<std::string, double> breads{{"wheat", 4.0}, {"white", 2.5}};
 * breakfast["toast"] = breads;
 * // Equivalent to:
 * breakfast["toast"]["wheat"] = 4.0
 * breakfast["toast"]["white"] = 2.5
 * ```
 *
 * ## Accessing elements
 *
 * ```
 * double val1 = breakfast["spam"].asDouble();
 * std::string val2 = breakfast["eggs"].asString();
 * vector_fp val3 = breakfast["beans"]["baked"].asVector<double>();
 *
 * std::map<std::string, double> = breakfast["toast"].asMap<double>();
 * ```
 *
 * ## Checking for elements
 *
 * ```
 * try {
 *     breakfast["waffle"].asDouble();
 * } except (std::exception& err) {
 *     // Exception will be thrown.
 *     // 'breakfast' will have an empty key named "waffle"
 * }
 *
 * try {
 *     breakfast.at("grits").asDouble();
 * } except (std::exception& err) {
 *     // Exception will be thrown and no new key will be added
 * }
 *
 * if (breakfast.hasKey("grits")) {
 *     // do something with this entry
 * }
 * ```
 *
 * ## Checking element types
 *
 * ```
 * if (breakfast["sausage"].is<vector<double>>()) {
 *     // access using asVector<double>
 * } else if (breakfast["sausage"].type() == typeid(vector<std::string>)) {
 *     // access using asVector<std::string>
 * }
 * ```
 */
class AnyMap
{
public:
    AnyMap() {};

    AnyValue& operator[](const std::string& key);

    AnyValue& at(const std::string& key);

    bool hasKey(const std::string& key) const;

private:
    std::unordered_map<std::string, AnyValue> m_data;
    friend class AnyValue;
};

// Definitions for templated functions which require the full declaration of
// class AnyMap.

template<class T>
AnyValue& AnyValue::operator=(const std::unordered_map<std::string, T> items) {
    m_value = AnyMap();
    AnyMap& dest = as<AnyMap>();
    for (const auto& item : items) {
        dest[item.first] = item.second;
    }
    return *this;
}

template<class T>
AnyValue& AnyValue::operator=(const std::map<std::string, T> items) {
    m_value = AnyMap();
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
        if (m_value.type() == typeid(void)) {
            m_value = std::move(AnyMap());
        }
        return boost::any_cast<AnyMap&>(m_value);
    } catch (boost::bad_any_cast&) {
        throw CanteraError("AnyValue::as",
            "value of key '{}' is a '{}',\nnot an 'AnyMap'.",
            m_key, demangle(m_value.type()));
    }
}

template<class T>
std::map<std::string, T> AnyValue::asMap()
{
    std::map<std::string, T> dest;
    for (const auto& item : as<AnyMap>().m_data) {
        try {
            dest[item.first] = boost::any_cast<T>(item.second.m_value);
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
