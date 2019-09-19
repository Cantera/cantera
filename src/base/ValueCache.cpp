/**
 *  @file ValueCache.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ValueCache.h"
#include <mutex>

namespace
{
std::mutex id_mutex;
}

namespace Cantera
{

int ValueCache::m_last_id = 0;

int ValueCache::getId()
{
    std::unique_lock<std::mutex> lock(id_mutex);
    return ++m_last_id;
}

void ValueCache::clear()
{
    m_scalarCache.clear();
    m_arrayCache.clear();
}

}
