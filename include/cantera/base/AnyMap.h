//! @file AnyMap.h

#ifndef CT_ANYMAP_H
#define CT_ANYMAP_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/Units.h"
#include "cantera/base/ctexceptions.h"

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>

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
    const AnyValue& operator[](const std::string& key) const;

    bool hasKey(const std::string& key) const;

    // The value knows the name of its corresponding key in order to provide
    // comprehensible error messages.
    void setKey(const std::string& key);

    template<class T>
    const T& as() const;

    template<class T>
    T& as();

    const std::type_info& type() const;
    std::string type_str() const;

    template<class T>
    bool is() const;

    explicit AnyValue(const std::string& value);
    explicit AnyValue(const char* value);
    AnyValue& operator=(const std::string& value);
    AnyValue& operator=(const char* value);
    const std::string& asString() const;

    explicit AnyValue(double value);
    AnyValue& operator=(double value);
    double& asDouble();
    const double& asDouble() const;

    explicit AnyValue(bool value);
    AnyValue& operator=(bool value);
    bool& asBool();
    const bool& asBool() const;

    explicit AnyValue(long int value);
    AnyValue& operator=(long int value);
    AnyValue& operator=(int value);
    long int& asInt();
    const long int& asInt() const;

    template<class T>
    AnyValue& operator=(const std::vector<T>& value);
    template<class T>
    const std::vector<T>& asVector(size_t nMin=npos, size_t nMax=npos) const;
    template<class T>
    std::vector<T>& asVector(size_t nMin=npos, size_t nMax=npos);

    AnyValue& operator=(const AnyMap& value);
    AnyValue& operator=(AnyMap&& value);

    template<class T>
    AnyValue& operator=(const std::unordered_map<std::string, T> items);

    template<class T>
    AnyValue& operator=(const std::map<std::string, T> items);

    template<class T>
    std::map<std::string, T> asMap() const;

    //! Access a vector<AnyMap> as a mapping using the value of `name` from each
    //! item as the key in the new mapping.
    /*!
     * For example, for the list:
     * ```
     * [{name: O2, weight: 32}, {name: CH4, weight: 16}]
     * ```
     * calling `asMap("name")` will create a map with keys ``O2`` and ``CH4``.
     */
    std::unordered_map<std::string, const AnyMap*> asMap(const std::string& name) const;
    std::unordered_map<std::string, AnyMap*> asMap(const std::string& name);

    //! @see AnyMap::applyUnits
    void applyUnits(const UnitSystem& units);

private:
    std::string demangle(const std::type_info& type) const;

    template<class T>
    void checkSize(const std::vector<T>& v, size_t nMin, size_t nMax) const;

    std::string m_key;
    std::unique_ptr<boost::any> m_value;
    static std::map<std::string, std::string> s_typenames;
};

// Implicit conversion to vector<AnyValue>
template<>
const std::vector<AnyValue>& AnyValue::asVector<AnyValue>(size_t nMin, size_t nMax) const;

template<>
std::vector<AnyValue>& AnyValue::asVector<AnyValue>(size_t nMin, size_t nMax);

// Implicit conversion of long int to double if accessed as a vector<double>
template<>
const std::vector<double>& AnyValue::asVector<double>(size_t nMin, size_t nMax) const;

template<>
std::vector<double>& AnyValue::asVector<double>(size_t nMin, size_t nMax);

// Implicit conversion of long int to double if accessed as a vector<vector<double>>
template<>
const std::vector<vector_fp>& AnyValue::asVector<vector_fp>(size_t nMin, size_t nMax) const;

template<>
std::vector<vector_fp>& AnyValue::asVector<vector_fp>(size_t nMin, size_t nMax);


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
 *     // 'breakfast' will have an empty key named "waffle" unless `breakfast`
 *     // is a `const AnyMap`.
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

    //! Create an AnyMap from a YAML file.
    /*!
     *  Searches the directory containing the optionally-specified parent file
     *  first, followed by the current working directory and the Cantera include
     *  path.
     */
    static AnyMap fromYamlFile(const std::string& name,
                               const std::string& parent_name="");

    //! Create an AnyMap from a string containing a YAML document
    static AnyMap fromYamlString(const std::string& yaml);

    AnyValue& operator[](const std::string& key);
    const AnyValue& operator[](const std::string& key) const;

    const AnyValue& at(const std::string& key) const;

    bool hasKey(const std::string& key) const;

    void erase(const std::string& key);

    //! Return a string listing the keys in this AnyMap, e.g. for use in error
    //! messages
    std::string keys_str() const;

    bool getBool(const std::string& key, bool default_) const;
    long int getInt(const std::string& key, long int default_) const;
    double getDouble(const std::string& key, double default_) const;
    const std::string& getString(const std::string& key,
                                 const std::string& default_) const;

    //! Convert the item stored by the given `key` to the units specified in
    //! `units`. If the stored value is a double, convert it using the default
    //! units. If the input is a string, treat this as a dimensioned value, e.g.
    //! '988 kg/m^3' and convert from the specified units.
    double convert(const std::string& key, const std::string& units) const;

    //! Convert the item stored by the given `key` to the units specified in
    //! `units`. If the stored value is a double, convert it using the default
    //! units. If the input is a string, treat this as a dimensioned value, e.g.
    //! '988 kg/m^3' and convert from the specified units. If the key is
    //! missing, the `default_` value is returned.
    double convert(const std::string& key, const std::string& units,
                   double default_) const;

    //! Convert a vector of dimensional values
    /*!
     * For each item in the vector, if the stored value is a double, convert it
     * using the default units. If the value is a string, treat it as a
     * dimensioned value, e.g. '988 kg/m^3', and convert from the specified
     * units.
     *
     * @param key    Location of the vector in this AnyMap
     * @param units  Units to convert to
     * @param nMin   Minimum allowed length of the vector. If #nMax is not
     *     specified, this is also taken to be the maximum length. An exception
     *     is thrown if this condition is not met.
     * @param nMax   Maximum allowed length of the vector. An exception is
     *     thrown if this condition is not met.
     */
    vector_fp convertVector(const std::string& key, const std::string& units,
                            size_t nMin=npos, size_t nMax=npos) const;

    //! Convert the item stored by the given `key` to the units specified in
    //! `units`. If the stored value is a double, convert it using the default
    //! units. If the input is a string, treat this as a dimensioned value, e.g.
    //! '2.7e4 J/kmol' and convert from the specified units.
    double convertMolarEnergy(const std::string& key,
                              const std::string& units) const;

    //! Convert the item stored by the given `key` to the units specified in
    //! `units`. If the stored value is a double, convert it using the default
    //! units. If the stored value is a string, treat it as a dimensioned value,
    //! e.g. '2.7e4 J/kmol' and convert from the specified units. If the key is
    //! missing, the `default_` value is returned.
    double convertMolarEnergy(const std::string& key, const std::string& units,
                              double default_) const;

    // Define begin() and end() to allow use with range-based for loops
    using const_iterator = std::unordered_map<std::string, AnyValue>::const_iterator;
    const_iterator begin() const {
        return m_data.begin();
    }

    const_iterator end() const {
        return m_data.end();
    }

    size_t size() {
        return m_data.size();
    };

    //! Return the default units that should be used to convert stored values
    const UnitSystem& units() const { return m_units; }

    //! Use the supplied UnitSystem to set the default units, and recursively
    //! process overrides from nodes named `units`.
    /*!
     * If a `units` node is present in a map that contains other keys, the
     * specified units are taken to be the defaults for that map. If the map
     * contains only a `units` node, and is the first item in a list of maps,
     * then the specified units are taken to be the defaults for all the maps in
     * the list.
     *
     * After being processed, the `units` nodes are removed, so this function
     * should be called only once, on the root AnyMap. This function is called
     * automatically by the fromYamlFile() and fromYamlString() constructors.
     */
    void applyUnits(const UnitSystem& units);

private:
    template <class T>
    const T& get(const std::string& key, const T& default_,
                 std::function<const T&(const AnyValue*)> getter) const;

    std::unordered_map<std::string, AnyValue> m_data;
    UnitSystem m_units;

    //! Cache for previously-parsed input (YAML) files. The key is the full path
    //! to the file, and the second element of the value is the last-modified
    //! time for the file, which is used to enable change detection.
    static std::unordered_map<std::string, std::pair<AnyMap, int>> s_cache;

    friend class AnyValue;
};

// Define begin() and end() to allow use with range-based for loops
AnyMap::const_iterator begin(const AnyValue& v);
AnyMap::const_iterator end(const AnyValue& v);

}

#ifndef CANTERA_API_NO_BOOST
#include "cantera/base/AnyMap.inl.h"
#endif

#endif
