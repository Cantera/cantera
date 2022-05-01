//! @file AnyMap.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/AnyMap.h"
#include "application.h"
#include "cantera/base/yaml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"
#include "cantera/base/utilities.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <mutex>
#include <unordered_set>

namespace ba = boost::algorithm;
using std::vector;
using std::string;

namespace { // helper functions

std::mutex yaml_cache_mutex;
std::mutex yaml_field_order_mutex;
using namespace Cantera;

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
    int numDigit = 0;
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
            numDigit++;
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
            if (numExp > 1 || numDigit == 0 || i == str.size() - 1) {
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

long int getPrecision(const Cantera::AnyValue& precisionSource)
{
    long int precision = 15;
    auto& userPrecision = precisionSource.getMetadata("precision");
    if (userPrecision.is<long int>()) {
        precision = userPrecision.asInt();
    }
    return precision;
}

string formatDouble(double x, long int precision)
{
    // This function ensures that trailing zeros resulting from round-off error
    // are removed. Values are only rounded if at least three digits are removed,
    // or the displayed value has multiple trailing zeros.
    if (x == 0.0) {
        return "0.0";
    }

    // Build string with full precision
    bool useExp = std::abs(x) < 1e-2 || std::abs(x) >= 1e4;
    int log10x;
    size_t last;
    string s0;
    if (useExp) {
        s0 = fmt::format(fmt::format("{:.{}e}", x, precision));
        // last digit of significand
        last = s0.size() - 5;
        if (s0[last + 1] == 'e') {
            // pass - most values use four letter exponent (examples: e+05, e-03)
        } else if (s0[last] == 'e') {
            last--; // exponents larger than e+99 or smaller than e-99 (example: e+100)
        } else {
            last = s0.find('e') - 1; // backstop; slower, but will always work
        }
    } else {
        log10x = static_cast<int>(std::floor(std::log10(std::abs(x))));
        s0 = fmt::format("{:.{}f}", x, precision - log10x);
        last = s0.size() - 1; // last digit
    }
    if (s0[last - 2] == '0' && s0[last - 1] == '0' && s0[last] < '5') {
        // Value ending in '00x' and should be rounded down
    } else if (s0[last - 2] == '9' && s0[last - 1] == '9' && s0[last] > '4') {
        // Value ending in '99y' and should be rounded up
    } else if (s0[last - 1] == '0' && s0[last] == '0') {
        // Remove trailing zeros
    } else {
        // Value should not be rounded / do not round last digit
        return s0;
    }

    // Remove trailing zeros
    string s1;
    if (s0[last - 1] == '0') {
        s1 = s0; // Recycle original string
    } else if (useExp) {
        s1 = fmt::format(fmt::format("{:.{}e}", x, precision - 2));
    } else {
        s1 = fmt::format("{:.{}f}", x, precision - log10x - 2);
    }
    size_t digit = last - 2;
    while (s1[digit] == '0' && s1[digit - 1] != '.') {
        digit--;
    }

    // Assemble rounded value and return
    if (useExp) {
        size_t eloc = s1.find('e');
        s0 = string(s1.begin() + eloc, s1.end());
    }
    s1 = string(s1.begin(), s1.begin() + digit + 1);
    if (useExp) {
        return s1 + s0;
    }
    return s1;
}

struct Quantity
{
    AnyValue value;
    Units units;
    bool isActivationEnergy;
    AnyValue::unitConverter converter;

    bool operator==(const Quantity& other) const {
    return value == other.value && units == other.units
            && isActivationEnergy == other.isActivationEnergy;
    }
};

Cantera::AnyValue Empty;

} // end anonymous namespace

namespace YAML { // YAML converters

using namespace Cantera;
static const int max_line_length = 87;

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
            AnyValue& value = target.createForYaml(key, loc.line, loc.column);
            if (child.second.IsMap()) {
                value = child.second.as<AnyMap>();
            } else {
                value = child.second.as<AnyValue>();
                value.setKey(key);
            }
        }
        return true;
    }
};

YAML::Emitter& operator<<(YAML::Emitter& out, const AnyMap& rhs)
{
    bool flow = rhs.getBool("__flow__", false);
    if (flow) {
        out << YAML::Flow;
        out << YAML::BeginMap;
        size_t width = 15;
        for (const auto& item : rhs.ordered()) {
            const auto& name = item.first;
            const auto& value = item.second;
            string valueStr;
            bool foundType = true;
            if (value.is<double>()) {
                valueStr = formatDouble(value.asDouble(), getPrecision(value));
            } else if (value.is<string>()) {
                valueStr = value.asString();
            } else if (value.is<long int>()) {
                valueStr = fmt::format("{}", value.asInt());
            } else if (value.is<bool>()) {
                valueStr = fmt::format("{}", value.asBool());
            } else {
                foundType = false;
            }

            if (foundType) {
                // Check if this item will fit on the current line, including spaces
                // for delimiters and whitespace. If not, wrap to the next line.
                if (width + name.size() + valueStr.size() + 4 > max_line_length) {
                    out << YAML::Newline;
                    width = 15;
                }
                out << name;
                out << valueStr;
                width += name.size() + valueStr.size() + 4;
            } else {
                // Put items of an unknown (compound) type on a line alone
                out << YAML::Newline;
                out << name;
                out << value;
                width = 99; // Force newline after this item as well
            }
        }
    } else {
        out << YAML::BeginMap;
        for (const auto& item : rhs.ordered()) {
            out << item.first;
            out << item.second;
        }
    }
    out << YAML::EndMap;
    return out;
}

//! Write YAML strings spanning multiple lines if input includes endline '\n'
void emitString(YAML::Emitter& out, const string& str0) {
    size_t endline = str0.rfind('\n');
    if (endline == std::string::npos) {
        out << str0;
        return;
    }

    // Remove trailing line break
    string str1 = str0;
    if (endline == str1.size() - 1) {
        str1.erase(endline, 1);
        endline = str1.rfind('\n');
    }

    // Deblank lines (remove whitespace surrounding line breaks)
    while (endline != std::string::npos) {
        size_t len = 1;
        while (str1[endline + len] == ' ') {
            len++; // account for whitespace after line break
        }
        while (str1[endline - 1] == ' ') {
            len++; // account for whitespace before line break
            endline--;
        }
        if (len > 1) {
            // remove surrounding whitespace
            str1.replace(endline, len, "\n");
        }
        endline = str1.rfind('\n', endline - 1);
    }
    out << YAML::Literal << str1;
}

//! Write a vector in YAML "flow" style, wrapping lines to avoid exceeding the
//! preferred maximum line length (set by `max_line_length`). Specialized for
//! `vector<double>` to be able to use the custom `formatDouble` function with
//! a given precision.
void emitFlowVector(YAML::Emitter& out, const vector<double>& v, long int precision)
{
    out << YAML::Flow;
    out << YAML::BeginSeq;
    size_t width = 15; // wild guess, but no better value is available
    for (auto& x : v) {
        string xstr = formatDouble(x, precision);
        // Wrap to the next line if this item would exceed the target line length
        if (width + xstr.size() > max_line_length) {
            out << YAML::Newline;
            width = 15;
        }
        out << xstr;
        width += xstr.size() + 2; // Update width including comma and space
    }
    out << YAML::EndSeq;
}

//! Write a vector in YAML "flow" style, wrapping lines to avoid exceeding the
//! preferred maximum line length (set by `max_line_length`).
template <typename T>
void emitFlowVector(YAML::Emitter& out, const vector<T>& v)
{
    out << YAML::Flow;
    out << YAML::BeginSeq;
    size_t width = 15; // wild guess, but no better value is available
    for (const auto& x : v) {
        string xstr = fmt::format("{}", x);
        // Wrap to the next line if this item would exceed the target line length
        if (width + xstr.size() > max_line_length) {
            out << YAML::Newline;
            width = 15;
        }
        out << xstr;
        width += xstr.size() + 2;
    }
    out << YAML::EndSeq;
}

YAML::Emitter& operator<<(YAML::Emitter& out, const AnyValue& rhs)
{
    if (rhs.isScalar()) {
        if (rhs.is<string>()) {
            emitString(out, rhs.asString());
        } else if (rhs.is<double>()) {
            out << formatDouble(rhs.asDouble(), getPrecision(rhs));
        } else if (rhs.is<long int>()) {
            out << rhs.asInt();
        } else if (rhs.is<bool>()) {
            out << rhs.asBool();
        } else {
            throw CanteraError("operator<<(YAML::Emitter&, AnyValue&)",
                "Don't know how to encode value of type '{}' with key '{}'",
                rhs.type_str(), rhs.m_key);
        }
    } else if (rhs.is<AnyMap>()) {
        out << rhs.as<AnyMap>();
    } else if (rhs.is<vector<AnyMap>>()) {
        out << rhs.asVector<AnyMap>();
    } else if (rhs.is<vector<double>>()) {
        emitFlowVector(out, rhs.asVector<double>(), getPrecision(rhs));
    } else if (rhs.is<vector<string>>()) {
        emitFlowVector(out, rhs.asVector<string>());
    } else if (rhs.is<vector<long int>>()) {
        emitFlowVector(out, rhs.asVector<long int>());
    } else if (rhs.is<vector<bool>>()) {
        emitFlowVector(out, rhs.asVector<bool>());
    } else if (rhs.is<vector<Cantera::AnyValue>>()) {
        out << rhs.asVector<Cantera::AnyValue>();
    } else if (rhs.is<vector<vector<double>>>()) {
        const auto& v = rhs.asVector<vector<double>>();
        long int precision = getPrecision(rhs);
        out << YAML::BeginSeq;
        for (const auto& u : v) {
            emitFlowVector(out, u, precision);
        }
        out << YAML::EndSeq;
    } else if (rhs.is<vector<vector<string>>>()) {
        const auto& v = rhs.asVector<vector<string>>();
        out << YAML::BeginSeq;
        for (const auto& u : v) {
            emitFlowVector(out, u);
        }
        out << YAML::EndSeq;
    } else if (rhs.is<vector<vector<long int>>>()) {
        const auto& v = rhs.asVector<vector<long int>>();
        out << YAML::BeginSeq;
        for (const auto& u : v) {
            emitFlowVector(out, u);
        }
        out << YAML::EndSeq;
    } else if (rhs.is<vector<vector<bool>>>()) {
        const auto& v = rhs.asVector<vector<bool>>();
        out << YAML::BeginSeq;
        for (const auto& u : v) {
            emitFlowVector(out, u);
        }
        out << YAML::EndSeq;
    } else {
        throw CanteraError("operator<<(YAML::Emitter&, AnyValue&)",
            "Don't know how to encode value of type '{}' with key '{}'",
            rhs.type_str(), rhs.m_key);
    }
    return out;
}


template<>
struct convert<Cantera::AnyValue> {
    static Node encode(const Cantera::AnyValue& rhs) {
        throw NotImplementedError("");
    }

    static bool decode(const Node& node, Cantera::AnyValue& target) {
        target.setLoc(node.Mark().line, node.Mark().column);
        if (node.IsScalar()) {
            // Scalar nodes are int, doubles, or strings
            std::string nodestr = node.as<std::string>();
            if (node.Tag() == "!") {
                // Prevent quoted strings from being implicitly converted to
                // numeric types. For example, the quoted YAML string '12345' should not
                // be interpreted as an integer
                target = nodestr;
            } else if (isInt(nodestr)) {
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
                } else if (subtypes == Type::String) {
                    target = node.as<std::vector<std::vector<std::string>>>();
                } else if (subtypes == Type::Bool) {
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
        } else if (node.IsNull()) {
            target = Empty;
            return true;
        }
        return false;
    }
};

}

namespace Cantera {

std::unordered_map<std::string, std::pair<AnyMap, int>> AnyMap::s_cache;

std::unordered_map<std::string, std::vector<std::string>> AnyMap::s_headFields;
std::unordered_map<std::string, std::vector<std::string>> AnyMap::s_tailFields;

// Methods of class AnyBase

AnyBase::AnyBase()
    : m_line(-1)
    , m_column(0)
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
    , m_equals(other.m_equals)
{
}

AnyValue& AnyValue::operator=(AnyValue const& other) {
    if (this == &other) {
        return *this;
    }
    AnyBase::operator=(other);
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
    m_equals = other.m_equals;
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

bool AnyValue::empty() const {
    return m_value->empty();
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

// Specialization for "Quantity"

void AnyValue::setQuantity(double value, const std::string& units, bool is_act_energy) {
    *m_value = Quantity{AnyValue(value), Units(units), is_act_energy};
    m_equals = eq_comparer<Quantity>;
}

void AnyValue::setQuantity(double value, const Units& units) {
    *m_value = Quantity{AnyValue(value), units, false};
    m_equals = eq_comparer<Quantity>;
}

void AnyValue::setQuantity(const vector_fp& values, const std::string& units) {
    AnyValue v;
    v = values;
    *m_value = Quantity{v, Units(units), false};
    m_equals = eq_comparer<Quantity>;
}

void AnyValue::setQuantity(const AnyValue& value, const unitConverter& converter)
{
    *m_value = Quantity{value, Units(0.0), false, converter};
    m_equals = eq_comparer<Quantity>;
}

template<>
bool AnyValue::is<vector<double>>() const
{
    if (m_value->type() == typeid(vector<double>)) {
        return true;
    } else if (m_value->type() == typeid(vector<AnyValue>)) {
        for (const auto& item : as<vector<AnyValue>>()) {
            if (!(item.is<double>()
                || (item.is<Quantity>() && item.as<Quantity>().value.is<double>())))
            {
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
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
            "Element is not a mapping or list of mappings.\n"
            "Looking for a mapping with key '{}' = '{}'", key, value);
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
            "Element is not a mapping or list of mappings.\n"
            "Looking for a mapping with key '{}' = '{}'", key, value);
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

std::pair<int, int> AnyValue::order() const
{
    return {m_line, m_column};
}

void AnyValue::applyUnits(shared_ptr<UnitSystem>& units)
{
    if (is<AnyMap>()) {
        AnyMap& m = as<AnyMap>();

        if (m.getBool("__unconvertible__", false)) {
            AnyMap delta = units->getDelta(UnitSystem());
            if (delta.hasKey("length") || delta.hasKey("quantity")
                || delta.hasKey("time"))
            {
                throw CanteraError("AnyValue::applyUnits", "AnyMap contains values"
                    " that cannot be converted to non-default unit systems\n(probably"
                    " reaction rates not associated with a Kinetics object)");
            }
        }
        // Units declaration applicable to this map
        m.applyUnits(units);
    } else if (is<std::vector<AnyMap>>()) {
        auto& list = as<std::vector<AnyMap>>();
        if (list.size() && list[0].hasKey("units") && list[0].size() == 1) {
            // First item in the list is a units declaration, which applies to
            // the items in the list
            auto deltaUnits = list[0]["units"];
            list[0].m_data.erase("units");
            for (auto& item : list) {
                if (item.hasKey("units")) {
                    if (item.size() == 1) {
                        // Any additional units declarations are errors
                        throw InputFileError("AnyValue::applyUnits", item,
                            "Found units entry as not the first item in a list.");
                    } else {
                        // Merge with a child units declaration
                        auto& childUnits = item["units"].as<AnyMap>();
                        for (auto& jtem : deltaUnits) {
                            if (!childUnits.hasKey(jtem.first)) {
                                childUnits[jtem.first] = jtem.second;
                            }
                        }
                    }
                } else if (item.hasKey("__units__")) {
                    // Merge with a child units declaration
                    auto& childUnits = item["__units__"].as<AnyMap>();
                    for (auto& jtem : deltaUnits) {
                        if (!childUnits.hasKey(jtem.first)) {
                            childUnits[jtem.first] = jtem.second;
                        }
                    }
                } else {
                    item["__units__"] = deltaUnits;
                }
                item.applyUnits(units);
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
    } else if (is<vector<AnyValue>>()) {
        for (auto& v : as<vector<AnyValue>>()) {
            v.applyUnits(units);
        }
    } else if (is<Quantity>()) {
        auto& Q = as<Quantity>();
        if (Q.converter) {
            Q.converter(Q.value, *units);
            m_equals = Q.value.m_equals;
            // Replace the value last since Q is a reference to m_value and won't be
            // valid after this
            m_value = std::move(Q.value.m_value);
        } else if (Q.value.is<double>()) {
            if (Q.isActivationEnergy) {
                *this = Q.value.as<double>() / units->convertActivationEnergyTo(1.0, Q.units);
            } else {
                *this = Q.value.as<double>() / units->convertTo(1.0, Q.units);
            }
        } else if (Q.value.is<vector_fp>()) {
            double factor = 1.0 / units->convertTo(1.0, Q.units);
            auto& old = Q.value.asVector<double>();
            vector_fp converted(old.size());
            scale(old.begin(), old.end(), converted.begin(), factor);
            *this = std::move(converted);
        } else {
            throw CanteraError("AnyValue::applyUnits", "Don't know how to "
                "convert Quantity with held type '{}' in key '{}'",
                Q.value.type_str(), m_key);
        }
    }
}

void AnyValue::setFlowStyle(bool flow)
{
    as<AnyMap>().setFlowStyle();
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

AnyMap::AnyMap()
    : m_units(new UnitSystem())
{
}

AnyValue& AnyMap::operator[](const std::string& key)
{
    const auto& iter = m_data.find(key);
    if (iter == m_data.end()) {
        // Create a new key to return
        // NOTE: 'insert' can be replaced with 'emplace' after support for
        // G++ 4.7 is dropped.
        AnyValue& value = m_data.insert({key, AnyValue()}).first->second;
        value.setKey(key);
        if (m_metadata) {
            value.propagateMetadata(m_metadata);
        }

        // A pseudo-location used to set the ordering when outputting to
        // YAML so nodes added this way will come before nodes from YAML,
        // with insertion order preserved.
        value.setLoc(-1, m_column);
        m_column += 10;

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

AnyValue& AnyMap::createForYaml(const std::string& key, int line, int column)
{
    // NOTE: 'insert' can be replaced with 'emplace' after support for
    // G++ 4.7 is dropped.
    AnyValue& value = m_data.insert({key, AnyValue()}).first->second;
    value.setKey(key);
    if (m_metadata) {
        value.propagateMetadata(m_metadata);
    }

    value.setLoc(line, column);
    return value;
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

bool AnyMap::empty() const
{
    return m_data.size() == 0;
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

void AnyMap::update(const AnyMap& other, bool keepExisting)
{
    for (const auto& item : other) {
        if (!keepExisting || !hasKey(item.first)) {
            (*this)[item.first] = item.second;
        }
    }
}

std::string AnyMap::keys_str() const
{
    fmt::memory_buffer b;
    auto iter = this->begin();
    if (iter != this->end()) {
        fmt_append(b, "{}", iter->first);
        ++iter;
    }
    while (iter != this->end()) {
        fmt_append(b, ", {}", iter->first);
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

void AnyMap::copyMetadata(const AnyMap& other)
{
    m_line = other.m_line;
    m_column = other.m_column;
    if (!other.m_metadata) {
        return;
    }

    if (m_metadata) {
        // Fork the metadata tree at this point to avoid affecting parent nodes
        m_metadata = make_shared<AnyMap>(*m_metadata);
    } else {
        m_metadata = make_shared<AnyMap>();
    }

    for (const auto& item : *other.m_metadata) {
        (*m_metadata)[item.first] = item.second;
    }

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


AnyMap::OrderedProxy::OrderedProxy(const AnyMap& data)
    : m_data(&data)
{
    // Units always come first
    if (m_data->hasKey("__units__") && m_data->at("__units__").as<AnyMap>().size()) {
        m_units.reset(new std::pair<const string, AnyValue>{"units", m_data->at("__units__")});
        m_units->second.setFlowStyle();
        m_ordered.emplace_back(std::pair<int, int>{-2, 0}, m_units.get());
    }

    int head = 0; // sort key of the first programmatically-added item
    int tail = 0; // sort key of the last programmatically-added item
    for (auto& item : *m_data) {
        const auto& order = item.second.order();
        if (order.first == -1) { // Item is not from an input file
            head = std::min(head, order.second);
            tail = std::max(tail, order.second);
        }
        m_ordered.emplace_back(order, &item);
    }
    std::sort(m_ordered.begin(), m_ordered.end());

    // Adjust sort keys for items that should moved to the beginning or end of
    // the list
    if (m_data->hasKey("__type__")) {
        bool order_changed = false;
        const auto& itemType = m_data->at("__type__").asString();
        std::unique_lock<std::mutex> lock(yaml_field_order_mutex);
        if (AnyMap::s_headFields.count(itemType)) {
            for (const auto& key : AnyMap::s_headFields[itemType]) {
                for (auto& item : m_ordered) {
                    if (item.first.first >= 0) {
                        // This and following items come from an input file and
                        // should not be re-ordered
                        break;
                    }
                    if (item.second->first == key) {
                        item.first.second = --head;
                        order_changed = true;
                    }
                }
            }
        }
        if (AnyMap::s_tailFields.count(itemType)) {
            for (const auto& key : AnyMap::s_tailFields[itemType]) {
                for (auto& item : m_ordered) {
                    if (item.first.first >= 0) {
                        // This and following items come from an input file and
                        // should not be re-ordered
                        break;
                    }
                    if (item.second->first == key) {
                        item.first.second = ++tail;
                        order_changed = true;
                    }
                }
            }
        }

        if (order_changed) {
            std::sort(m_ordered.begin(), m_ordered.end());
        }
    }
}

AnyMap::OrderedIterator AnyMap::OrderedProxy::begin() const
{
    return OrderedIterator(m_ordered.begin(), m_ordered.end());
}

AnyMap::OrderedIterator AnyMap::OrderedProxy::end() const
{
    return OrderedIterator(m_ordered.end(), m_ordered.end());
}

AnyMap::OrderedIterator::OrderedIterator(
    const AnyMap::OrderedProxy::OrderVector::const_iterator& start,
    const AnyMap::OrderedProxy::OrderVector::const_iterator& stop)
{
    m_iter = start;
    m_stop = stop;
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

void AnyMap::applyUnits()
{
    applyUnits(m_units);
}

void AnyMap::applyUnits(shared_ptr<UnitSystem>& units) {
    if (hasKey("units")) {
        m_data["__units__"] = std::move(m_data["units"]);
        m_data.erase("units");
    }
    if (hasKey("__units__")) {
        m_units.reset(new UnitSystem(*units));
        m_units->setDefaults(m_data["__units__"].asMap<std::string>());
    } else {
        m_units = units;
    }
    for (auto& item : m_data) {
        item.second.applyUnits(m_units);
    }
}

void AnyMap::setUnits(const UnitSystem& units)
{
    if (hasKey("__units__")) {
        for (const auto& item : units.getDelta(*m_units)) {
            m_data["__units__"][item.first] = item.second;
        }
    } else {
        m_data["__units__"] = units.getDelta(*m_units);
    }
    m_units.reset(new UnitSystem(units));
}

void AnyMap::setFlowStyle(bool flow) {
    (*this)["__flow__"] = flow;
}

bool AnyMap::addOrderingRules(const string& objectType,
                             const vector<vector<string>>& specs)
{
    std::unique_lock<std::mutex> lock(yaml_field_order_mutex);
    for (const auto& spec : specs) {
        if (spec.at(0) == "head") {
            s_headFields[objectType].push_back(spec.at(1));
        } else if (spec.at(0) == "tail") {
            s_tailFields[objectType].push_back(spec.at(1));
        } else {
            throw CanteraError("AnyMap::addOrderingRules",
                "Unknown ordering rule '{}'", spec.at(0));
        }
    }
    return true;
}

void AnyMap::clearCachedFile(const std::string& filename)
{
    std::string fullName = findInputFile(filename);
    if (s_cache.count(fullName)) {
        s_cache.erase(fullName);
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
    amap.applyUnits();
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
        cache_item.first.applyUnits();
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

std::string AnyMap::toYamlString() const
{
    YAML::Emitter out;
    const_cast<AnyMap*>(this)->applyUnits();
    out << *this;
    out << YAML::Newline;
    return out.c_str();
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

    fmt_append(b, "|  Line |\n");
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
            fmt_append(b, "> {: 5d} > {}\n", i+1, line);
            fmt_append(b, "{:>{}}\n", "^", column + 11);
            lastShown = i;
        } else if ((lineno + 4 > i && lineno < i + 6) ||
                   (lineno2 + 4 > i && lineno2 < i + 6)) {
            if (lastShown >= 0 && i - lastShown > 1) {
                fmt_append(b, "...\n");
            }
            fmt_append(b, "| {: 5d} | {}\n", i+1, line);
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
    fmt_append(b, "Error on line {} of {}:\n{}\n", lineno+1, filename, message);
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
        fmt_append(b, "Error on lines {} and {} of {}:\n",
                   std::min(line1, line2) + 1, std::max(line1, line2) + 1, filename1);
        fmt_append(b, "{}\n", message);
        formatInputFile(b, metadata1, filename1, line1, column1, line2, column2);
    } else {
        fmt_append(b, "Error on line {} of {} and line {} of {}:\n{}\n",
                   line1+1, filename1, line2+1, filename2, message);
        formatInputFile(b, metadata1, filename1, line1, column1);
        fmt_append(b, "\n");
        formatInputFile(b, metadata2, filename2, line2, column2);
    }

    return to_string(b);
}

void warn_deprecated(const std::string& source, const AnyBase& node,
                     const std::string& message)
{
    if (!node.m_metadata) {
        warn_deprecated(source, message);
        return;
    }

    std::string filename = node.m_metadata->getString("filename", "input string");
    fmt::memory_buffer b;
    fmt_append(b, message);
    fmt_append(b, "\n");
    fmt_append(b, "On line {} of {}:\n", node.m_line+1, filename);
    formatInputFile(b, node.m_metadata, filename, node.m_line, node.m_column);
    warn_deprecated(source, to_string(b));
}

}
