//! @file AnyMap.inl.h

#ifndef CT_ANYMAP_INL_H
#define CT_ANYMAP_INL_H

#include "cantera/base/AnyMap.h"

namespace Cantera
{
// re-declared to avoid needing to include global.h here
string demangle(const std::type_info& type);

// Definitions for AnyValue templated functions

template<class T>
const T &AnyValue::as() const {
    try {
        if (typeid(T) == typeid(double) && m_value.type() == typeid(long int)) {
            // Implicit conversion of long int to double
            const_cast<AnyValue*>(this)->m_value = static_cast<double>(as<long int>());
            m_equals = eq_comparer<double>;
        } else if (typeid(T) == typeid(string) && m_value.type() == typeid(double)) {
            // Implicit conversion of double to string
            const_cast<AnyValue*>(this)->m_value = fmt::format("{}", as<double>());
            m_equals = eq_comparer<string>;
        } else if (typeid(T) == typeid(string) && m_value.type() == typeid(long int)) {
            // Implicit conversion of long int to string
            const_cast<AnyValue*>(this)->m_value = fmt::format("{}", as<long int>());
            m_equals = eq_comparer<string>;
        } else if (typeid(T) == typeid(vector<double>)
                   && m_value.type() == typeid(vector<AnyValue>)) {
            // Implicit conversion of vector<AnyValue> to vector<double>
            auto& asAny = as<vector<AnyValue>>();
            vector<double> asDouble(asAny.size());
            for (size_t i = 0; i < asAny.size(); i++) {
                asDouble[i] = asAny[i].as<double>();
            }
            const_cast<AnyValue*>(this)->m_value = std::move(asDouble);
            m_equals = eq_comparer<vector<double>>;
        }
        return std::any_cast<const T&>(m_value);
    } catch (std::bad_any_cast&) {
        if (m_value.type() == typeid(void)) {
            // Values that have not been set are of type 'void'
            throw InputFileError("AnyValue::as", *this,
                "Key '{}' not found or contains no value", m_key);
        } else {
            if (m_key == "") {
                throw InputFileError("AnyValue::as", *this,
                    "Unable to convert '{}' to '{}'.",
                    demangle(m_value.type()), demangle(typeid(T)));
            }
            throw InputFileError("AnyValue::as", *this,
                "Key '{}' contains a '{}',\nnot a '{}'",
                m_key, demangle(m_value.type()), demangle(typeid(T)));
        }
    }
}

template<class T>
T &AnyValue::as() {
    // To avoid duplicating the code from the const version, call that version
    // and just remove the const specifier from the return value
    return const_cast<T&>(const_cast<const AnyValue*>(this)->as<T>());
}

template<class T>
bool AnyValue::is() const {
    return m_value.type() == typeid(T);
}

template<> bool AnyValue::is<vector<double>>() const;

template<class T>
bool AnyValue::isVector() const {
    return m_value.type() == typeid(vector<T>);
}

template<class T>
bool AnyValue::isMatrix(size_t cols) const {
    if (m_value.type() != typeid(vector<vector<T>>)) {
        // not a matrix
        return false;
    }
    auto& asMatrix = as<vector<vector<T>>>();
    if (!asMatrix.size()) {
        // empty matrix
        return true;
    }
    if (cols == npos) {
        cols = asMatrix[0].size();
    }
    for (const auto& row : asMatrix) {
        if (row.size() != cols) {
            return false;
        }
    }
    return true;
}

template<class T>
AnyValue &AnyValue::operator=(const vector<T> &value) {
    m_value = value;
    m_equals = eq_comparer<vector<T>>;
    return *this;
}

template<class T>
const vector<T> &AnyValue::asVector(size_t nMin, size_t nMax) const {
    const auto& v = as<vector<T>>();
    checkSize(v, nMin, nMax);
    return v;
}

template<class T>
vector<T> &AnyValue::asVector(size_t nMin, size_t nMax) {
    auto& v = as<vector<T>>();
    checkSize(v, nMin, nMax);
    return v;
}

template<class T>
AnyValue& AnyValue::operator=(const std::unordered_map<string, T> items) {
    m_value = AnyMap();
    m_equals = eq_comparer<AnyMap>;
    AnyMap& dest = as<AnyMap>();
    for (const auto& [key, value] : items) {
        dest[key] = value;
    }
    return *this;
}

template<class T>
AnyValue& AnyValue::operator=(const map<string, T> items) {
    m_value = AnyMap();
    m_equals = eq_comparer<AnyMap>;
    AnyMap& dest = as<AnyMap>();
    for (const auto& [key, value] : items) {
        dest[key] = value;
    }
    return *this;
}

template<>
inline AnyMap& AnyValue::as<AnyMap>() {
    try {
        // This is where nested AnyMaps are created when the syntax
        // m[key1][key2] is used.
        if (m_value.type() == typeid(void)) {
            m_value = AnyMap();
            m_equals = eq_comparer<AnyMap>;
        }
        return std::any_cast<AnyMap&>(m_value);
    } catch (std::bad_any_cast&) {
        throw InputFileError("AnyValue::as", *this,
            "value of key '{}' is a '{}',\nnot an 'AnyMap'.",
            m_key, demangle(m_value.type()));
    }
}

template<class T>
map<string, T> AnyValue::asMap() const
{
    map<string, T> dest;
    for (const auto& item : as<AnyMap>()) {
        dest[item.first] = item.second.as<T>();
    }
    return dest;
}

template<class T>
void AnyValue::checkSize(const vector<T>& v, size_t nMin, size_t nMax) const
{
    if (nMin != npos && nMax == npos && v.size() != nMin) {
        throw InputFileError("AnyValue::checkSize", *this,
            "Expected array '{}' to have length {}, but found "
            "an array of length {}.", m_key, nMin, v.size());
    } else if (nMin != npos && nMax != npos
               && (v.size() < nMin || v.size() > nMax)) {
        throw InputFileError("AnyValue::checkSize", *this,
            "Expected array '{}' to have from {} to {} elements, but found "
            "an array of length {}.", m_key, nMin, nMax, v.size());
    }
}

template<class T, class U>
bool AnyValue::vector_eq(const std::any& lhs, const std::any& rhs)
{
    const auto& lvec = std::any_cast<T>(lhs);
    const auto& rvec = std::any_cast<U>(rhs);
    if (lvec.size() != rvec.size()) {
        return false;
    } else {
        return std::equal(lvec.begin(), lvec.end(), rvec.begin());
    }
}

template<class T, class U>
bool AnyValue::vector2_eq(const std::any& lhs, const std::any& rhs)
{
    const auto& lvec = std::any_cast<vector<T>>(lhs);
    const auto& rvec = std::any_cast<vector<U>>(rhs);
    if (lvec.size() != rvec.size()) {
        return false;
    } else {
        for (size_t i = 0; i < lvec.size(); i++) {
            if (!std::equal(lvec[i].begin(), lvec[i].end(), rvec[i].begin())) {
                return false;
            }
        }
        return true;
    }
}

template<class T>
bool AnyValue::eq_comparer(const std::any& lhs, const std::any& rhs)
{
    using std::any_cast;
    typedef vector<double> vd;
    typedef vector<long int> vi;
    typedef vector<AnyValue> va;
    typedef vector<string> vs;

    auto& ltype = lhs.type();
    auto& rtype = rhs.type();
    AssertThrowMsg(ltype == typeid(T),
        "AnyValue::eq_comparer", "Compare function ({}) does not match held type ({})",
        demangle(ltype), demangle(typeid(T)));

    if (ltype == rtype) {
        return any_cast<T>(lhs) == any_cast<T>(rhs);
    } else if (ltype == typeid(double) && rtype == typeid(long int)) {
        return any_cast<double>(lhs) == any_cast<long int>(rhs);
    } else if (ltype == typeid(long int) && rtype == typeid(double)) {
        return any_cast<long int>(lhs) == any_cast<double>(rhs);

    } else if (ltype == typeid(vd) && rtype == typeid(vi)) {
        return vector_eq<vd, vi>(lhs, rhs);
    } else if (ltype == typeid(vi) && rtype == typeid(vd)) {
        return vector_eq<vi, vd>(lhs, rhs);

    } else if (ltype == typeid(va)) {
        if (rtype == typeid(vd)) {
            return vector_eq<va, vd>(lhs, rhs);
        } else if (rtype == typeid(vi)) {
            return vector_eq<va, vi>(lhs, rhs);
        } else if (rtype == typeid(vs)) {
            return vector_eq<va, vs>(lhs, rhs);
        }
    } else if (rtype == typeid(va)) {
        if (ltype == typeid(vd)) {
            return vector_eq<vd, va>(lhs, rhs);
        } else if (ltype == typeid(vi)) {
            return vector_eq<vi, va>(lhs, rhs);
        } else if (ltype == typeid(vs)) {
            return vector_eq<vs, va>(lhs, rhs);
        }
    } else if (ltype == typeid(vector<vd>) && rtype == typeid(vector<vi>)) {
        return vector2_eq<vd, vi>(lhs, rhs);
    } else if (ltype == typeid(vector<vi>) && rtype == typeid(vector<vd>)) {
        return vector2_eq<vd, vi>(lhs, rhs);
    }
    return false;
}

}
#endif
