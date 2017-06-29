//! @file AnyMap.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/AnyMap.h"

namespace Cantera
{

std::map<std::string, std::string> AnyValue::s_typenames = {
    {typeid(double).name(), "double"},
    {typeid(std::string).name(), "string"},
    {typeid(std::vector<double>).name(), "vector<double>"},
    {typeid(AnyMap).name(), "AnyMap"},
};

// Methods of class AnyValue

AnyValue& AnyValue::operator[](const std::string& key)
{
    return as<AnyMap>()[key];
}

bool AnyValue::hasKey(const std::string& key) const {
    return (is<AnyMap>() && as<AnyMap>().hasKey(key));
}

AnyValue& AnyValue::operator=(const AnyMap& value) {
    m_value = value;
    return *this;
}

AnyValue& AnyValue::operator=(AnyMap&& value) {
    m_value = std::move(value);
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
            value = std::move(AnyMap());
            value.setKey(head);
            return value.as<AnyMap>()[tail];
        } else {
            // Return an already existing key
            return iter->second.as<AnyMap>()[tail];
        }
    }
}

AnyValue& AnyMap::at(const std::string& key)
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

}
