/**
 *  @file ValueCache.cpp
 */

#include "cantera/base/ValueCache.h"
#include "cantera/base/ct_thread.h"

namespace
{
Cantera::mutex_t id_mutex;
}

namespace Cantera
{

int ValueCache::m_last_id = 0;

int ValueCache::getId()
{
    ScopedLock lock(id_mutex);
    return ++m_last_id;
}

void ValueCache::clear()
{
    m_scalarCache.clear();
    m_arrayCache.clear();
}

}
