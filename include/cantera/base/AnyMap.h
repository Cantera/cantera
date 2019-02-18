//! @file AnyMap.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

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
class InputFile;

//! A wrapper for a variable whose type is determined at runtime
/*!
 * Instances of AnyValue are used as values in an AnyMap. Values are converted
 * to a concrete type using the templated as() method or convenience methods
 * such as asString() and asDouble(). See AnyMap for usage examples.
 *
 * Elements are set using assignment, and the assignment operator has been
 * overloaded for specific types so that only those types are allowed to be
 * used in an AnyValue. The allowed types are:
 * - AnyMap
 * - `double`
 * - `long int`
 * - `bool`
 * - `std::string`
 * - `std::vector` of any of the above
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

    //! If this AnyValue is an AnyMap, return the value stored in `key`.
    AnyValue& operator[](const std::string& key);
    const AnyValue& operator[](const std::string& key) const;

    //! Returns `true` if this AnyValue is an AnyMap and that map contains
    //! a key with the given name.
    bool hasKey(const std::string& key) const;

    //! Set the name of the key storing this value in an AnyMap. Used for
    //! providing informative error messages in class InputFileError.
    void setKey(const std::string& key);

    //! For values which are derived from an input file, set the line and column
    //! of this value in that file. Used for providing context for some error
    //! messages.
    void setLoc(int line, int column);

    //! Set information about the file used to create this value. Recursively
    //! sets the information on any child elements.
    void setFile(shared_ptr<InputFile>& file);

    //! Get the value of this key as the specified type.
    template<class T>
    const T& as() const;

    template<class T>
    T& as();

    //! Returns the type of the held value.
    const std::type_info& type() const;

    //! Returns a string specifying the type of the held value.
    std::string type_str() const;

    //! Returns `true` if the held value is of the specified type.
    template<class T>
    bool is() const;

    //! Returns `true` if the held value is a scalar type (e.g. `double`, `long
    //! int`, `string`, or `bool`).
    bool isScalar() const;

    explicit AnyValue(const std::string& value);
    explicit AnyValue(const char* value);
    AnyValue& operator=(const std::string& value);
    AnyValue& operator=(const char* value);
    //! Return the held value, if it is a string
    const std::string& asString() const;

    explicit AnyValue(double value);
    AnyValue& operator=(double value);
    //! Return the held value as a `double`, if it is a `double` or a `long
    //! int`.
    double& asDouble();
    const double& asDouble() const;

    explicit AnyValue(bool value);
    AnyValue& operator=(bool value);
    //! Return the held value, if it is a `bool`.
    bool& asBool();
    const bool& asBool() const;

    explicit AnyValue(long int value);
    AnyValue& operator=(long int value);
    AnyValue& operator=(int value);
    //! Return the held value, if it is a `long int`.
    long int& asInt();
    const long int& asInt() const;

    template<class T>
    AnyValue& operator=(const std::vector<T>& value);
    //! Return the held value, if it is a vector of type `T`. If called with one
    //! argument, requires the vector to be of the specified size. If called
    //! with two arguments, requires the vector to be within the range specified
    //! by the two values, inclusive.
    template<class T>
    const std::vector<T>& asVector(size_t nMin=npos, size_t nMax=npos) const;
    template<class T>
    std::vector<T>& asVector(size_t nMin=npos, size_t nMax=npos);

    explicit AnyValue(const AnyMap& value);
    AnyValue& operator=(const AnyMap& value);
    AnyValue& operator=(AnyMap&& value);

    template<class T>
    AnyValue& operator=(const std::unordered_map<std::string, T> items);

    template<class T>
    AnyValue& operator=(const std::map<std::string, T> items);

    //! Return the held `AnyMap` as a `std::map` where all of the values have
    //! the specified type.
    template<class T>
    std::map<std::string, T> asMap() const;

    //! Access a `vector<AnyMap>` as a mapping using the value of `name` from
    //! each item as the key in the new mapping.
    /*!
     * For example, for the list:
     * ```
     * [{name: O2, weight: 32}, {name: CH4, weight: 16}]
     * ```
     * calling `asMap("name")` will create a map with keys ``O2`` and ``CH4``.
     */
    std::unordered_map<std::string, const AnyMap*> asMap(const std::string& name) const;
    std::unordered_map<std::string, AnyMap*> asMap(const std::string& name);

    //! For objects of type `vector<AnyMap>`, return the item where the given
    //! key has the specified value. If value is the empty string, returns the
    //! first item in the list.
    AnyMap& getMapWhere(const std::string& key, const std::string& value);

    //! @see AnyMap::applyUnits
    void applyUnits(const UnitSystem& units);

private:
    std::string demangle(const std::type_info& type) const;

    template<class T>
    void checkSize(const std::vector<T>& v, size_t nMin, size_t nMax) const;

    //! Line where this value occurs in the input file
    int m_line;

    //! Column where this value occurs in the input file
    int m_column;

    //! Information about the input file used to create this object
    shared_ptr<InputFile> m_file;

    //! Key of this value in a parent `AnyMap`
    std::string m_key;

    //! The held value
    std::unique_ptr<boost::any> m_value;

    //! Human-readable names for some common types, for use when
    //! `boost::demangle` is not available.
    static std::map<std::string, std::string> s_typenames;

    friend class InputFileError;
};

//! Implicit conversion to vector<AnyValue>
template<>
const std::vector<AnyValue>& AnyValue::asVector<AnyValue>(size_t nMin, size_t nMax) const;

template<>
std::vector<AnyValue>& AnyValue::asVector<AnyValue>(size_t nMin, size_t nMax);

//! Implicit conversion of long int to double if accessed as a vector<double>
template<>
const std::vector<double>& AnyValue::asVector<double>(size_t nMin, size_t nMax) const;

template<>
std::vector<double>& AnyValue::asVector<double>(size_t nMin, size_t nMax);

//! Implicit conversion of long int to double if accessed as a vector<vector<double>>
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

    //! Get the value of the item stored in `key`.
    AnyValue& operator[](const std::string& key);
    const AnyValue& operator[](const std::string& key) const;

    //! Get the value of the item stored in `key`. Raises an exception if the
    //! value does not exist.
    const AnyValue& at(const std::string& key) const;

    //! Returns `true` if the map contains an item named `key`.
    bool hasKey(const std::string& key) const;

    //! Erase the value held by `key`.
    void erase(const std::string& key);

    //! Return a string listing the keys in this AnyMap, e.g. for use in error
    //! messages
    std::string keys_str() const;

    //! For AnyMaps which are derived from an input file, set the line and
    //! column of this AnyMap in that file. Used for providing context for some
    //! error messages.
    void setLoc(int line, int column);

    //! Set information about the file used to create this AnyMap. Recursively
    //! sets the information on any child elements.
    void setFile(shared_ptr<InputFile>& file);

    //! Set the name of the file used to create this AnyMap. Recursively sets
    //! the information on any child elements.
    void setFileName(const std::string& filename);

    //! Set the contents of the file used to create this AnyMap. Used in the
    //! case where the AnyMap is created from an input string rather than a
    //! file. Recursively sets the information on any child elements.
    void setFileContents(const std::string& contents);

    //! If `key` exists, return it as a `bool`, otherwise return `default_`.
    bool getBool(const std::string& key, bool default_) const;

    //! If `key` exists, return it as a `long int`, otherwise return `default_`.
    long int getInt(const std::string& key, long int default_) const;

    //! If `key` exists, return it as a `double`, otherwise return `default_`.
    double getDouble(const std::string& key, double default_) const;

    //! If `key` exists, return it as a `string`, otherwise return `default_`.
    const std::string& getString(const std::string& key,
                                 const std::string& default_) const;

    //! Convert the item stored by the given `key` to the units specified in
    //! `units`. If the stored value is a double, convert it using the default
    //! units. If the input is a string, treat this as a dimensioned value, e.g.
    //! '988 kg/m^3' and convert from the specified units.
    double convert(const std::string& key, const std::string& units) const;
    double convert(const std::string& key, const Units& units) const;

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
     * @param nMin   Minimum allowed length of the vector. If `nMax` is not
     *     specified, this is also taken to be the maximum length. An exception
     *     is thrown if this condition is not met.
     * @param nMax   Maximum allowed length of the vector. An exception is
     *     thrown if this condition is not met.
     */
    vector_fp convertVector(const std::string& key, const std::string& units,
                            size_t nMin=npos, size_t nMax=npos) const;

    using const_iterator = std::unordered_map<std::string, AnyValue>::const_iterator;

    //! Defined to allow use with range-based for loops
    const_iterator begin() const {
        return m_data.begin();
    }

    //! Defined to allow use with range-based for loops
    const_iterator end() const {
        return m_data.end();
    }

    //! Returns the number of elements in this map
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
    //! The stored data
    std::unordered_map<std::string, AnyValue> m_data;

    //! The default units that are used to convert stored values
    UnitSystem m_units;

    //! Starting line for this map in the input file
    int m_line;

    //! Starting column for this map in the input file
    int m_column;

    //! Information about the file used to create this map
    shared_ptr<InputFile> m_file;

    //! Cache for previously-parsed input (YAML) files. The key is the full path
    //! to the file, and the second element of the value is the last-modified
    //! time for the file, which is used to enable change detection.
    static std::unordered_map<std::string, std::pair<AnyMap, int>> s_cache;

    friend class AnyValue;
    friend class InputFileError;
};

// Define begin() and end() to allow use with range-based for loops
AnyMap::const_iterator begin(const AnyValue& v);
AnyMap::const_iterator end(const AnyValue& v);

//! Error thrown for problems processing information contained in an AnyMap or
//! AnyValue.
/*!
 *  This class uses the file, line, and column information stored in an AnyMap
 *  or AnyValue to provide an error message including context lines for the
 *  original user input.
 */
class InputFileError : public CanteraError
{
public:
    //! Indicate an error occurring in `procedure` while using information from
    //! `node`. The `message` and `args` are processed as in the CanteraError
    //! class.
    template <typename... Args>
    InputFileError(const std::string& procedure, const AnyValue& node,
                   const std::string& message, const Args&... args)
        : CanteraError(
            procedure,
            formatError(fmt::format(message, args...),
                        node.m_line, node.m_column, node.m_file))
        {
        }

    //! Indicate an error occurring in `procedure` while using information from
    //! `node`. The `message` and `args` are processed as in the CanteraError
    //! class.
    template <typename... Args>
    InputFileError(const std::string& procedure, const AnyMap& node,
                   const std::string& message, const Args&... args)
        : CanteraError(
            procedure,
            formatError(fmt::format(message, args...),
                        node.m_line, node.m_column, node.m_file))
        {
        }

    virtual std::string getClass() const {
        return "InputFileError";
    }
protected:
    static std::string formatError(const std::string& message,
                                   int line, int column,
                                   const shared_ptr<InputFile>& file);
};

}

#ifndef CANTERA_API_NO_BOOST
#include "cantera/base/AnyMap.inl.h"
#endif

#endif
