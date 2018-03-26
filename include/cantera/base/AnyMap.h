//! @file AnyMap.h

#ifndef CT_ANYMAP_H
#define CT_ANYMAP_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

namespace boost
{
class any;
}

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
    AnyValue();
    ~AnyValue();
    AnyValue(AnyValue const& other);
    AnyValue(AnyValue&& other);
    AnyValue& operator=(AnyValue const& other);
    AnyValue& operator=(AnyValue&& other);

    AnyValue& operator[](const std::string& key);

    bool hasKey(const std::string& key) const;

    // The value knows the name of its corresponding key in order to provide
    // comprehensible error messages.
    void setKey(const std::string& key);

    template<class T>
    const T& as() const;

    template<class T>
    T& as();

    const std::type_info& type();

    template<class T>
    bool is() const;

    AnyValue& operator=(const std::string& value);
    AnyValue& operator=(const char* value);
    const std::string& asString() const;

    AnyValue& operator=(double value);
    double asDouble() const;

    AnyValue& operator=(bool value);
    bool asBool() const;

    AnyValue& operator=(long int value);
    AnyValue& operator=(int value);
    long int asInt() const;

    template<class T>
    AnyValue& operator=(const std::vector<T>& value);
    template<class T>
    const std::vector<T>& asVector() const;
    template<class T>
    std::vector<T>& asVector();

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
    std::unique_ptr<boost::any> m_value;
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

}

#ifndef CANTERA_API_NO_BOOST
#include "cantera/base/AnyMap.inl.h"
#endif

#endif
