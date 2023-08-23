/**
 *  @file ValueCache.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_VALUECACHE_H
#define CT_VALUECACHE_H

#include "ct_defs.h"
#include <limits>

namespace Cantera
{

/**
 * A cached property value and the state at which it was evaluated.
 *
 * This struct stores the value of some property evaluated at a particular
 * thermodynamic state. The #value can be either a real scalar or an array,
 * depending on the template parameter `T`. The exact meaning of #state1,
 * #state2, and #stateNum is determined by the function using the cached value,
 * which can check any combination of these variables before deciding whether
 * to recompute the cached values.
 *
 * References to CachedValue objects are returned by the "get" methods of
 * ValueCache, for example ValueCache::getScalar. Functions accessing cached values
 * should use the typedefs CachedScalar and CachedArray. See ValueCache for
 * details on how these classes should be used together.
 */
template <class T>
struct CachedValue {
    CachedValue() = default;

    //! Check whether the currently cached value is valid based on
    //! a single state variable. If it is not valid it updates the stored
    //! state to the new state in addition to returning false.
    bool validate(double state1New) {
      if(state1 == state1New) {
        return true;
      } else {
        state1 = state1New;
      }
      return false;
    }

    //! Check whether the currently cached value is valid based on
    //! state1 and state2. If it is not valid it updates the stored
    //! state to the new state in addition to returning false.
    bool validate(double state1New, double state2New) {
      if(state1 == state1New && state2 == state2New) {
        return true;
      } else {
        state1 = state1New;
        state2 = state2New;
      }
      return false;
    }

    //! Check whether the currently cached value is valid based on
    //! state1 and stateNum. If it is not valid it updates the stored
    //! state to the new state in addition to returning false.
    bool validate(double state1New, int stateNumNew) {
      if(state1 == state1New && stateNum == stateNumNew) {
        return true;
      } else {
        state1 = state1New;
        stateNum = stateNumNew;
      }
      return false;
    }

    //! Check whether the currently cached value is valid based on
    //! stateNum. If it is not valid it updates the stored
    //! state to the new state in addition to returning false.
    bool validate(int stateNumNew) {
      if(stateNum == stateNumNew) {
        return true;
      } else {
        stateNum = stateNumNew;
      }
      return false;
    }

    //! Check whether the currently cached value is valid based on
    //! state1, state2, and stateNum. If it is not valid it updates the stored
    //! state to the new state in addition to returning false.
    bool validate(double state1New, double state2New, int stateNumNew) {
      if(state1 == state1New && state2 == state2New && stateNum == stateNumNew) {
        return true;
      } else {
        state1 = state1New;
        state2 = state2New;
        stateNum = stateNumNew;
      }
      return false;
    }

    //! Value of the first state variable for the state at which #value was
    //! evaluated, for example temperature.
    double state1 = std::numeric_limits<double>::quiet_NaN();

    //! Value of the second state variable for the state at which #value was
    //! evaluated, for example density or pressure.
    double state2 = std::numeric_limits<double>::quiet_NaN();

    //! A surrogate for the composition. For cached properties of Phase,
    //! this should be set to Phase::stateMFNumber()
    int stateNum = std::numeric_limits<int>::min();

    //! The value of the cached property
    T value = T();
};

typedef CachedValue<double>& CachedScalar;
typedef CachedValue<vector<double>>& CachedArray;

/**
 * Storage for cached values.
 *
 * Stores cached values of properties evaluated at a particular thermodynamic
 * state. A class that needs cached values can have a ValueCache as a
 * member variable.
 *
 * Each method in the class that implements caching behavior needs a unique id
 * for its cached value. This id should be obtained by using the getId()
 * function to initialize a static variable within the function.
 *
 * For cases where the property is a scalar or vector, the cached value can be
 * stored in the CachedValue object. If the data type of the cached value is
 * more complex, then it can be stored in the calling class, and the value
 * attribute of the CachedScalar object can be ignored.
 *
 * An example use of class ValueCache:
 * @code
 * class Example {
 *     ValueCache m_cache;
 *     double get_property(double T, double P) {
 *         const static int cacheId = m_cache.getId();
 *         CachedScalar cached = m_cache.getScalar(cacheId);
 *         if (T != cached.state1 || P != cached.state2) {
 *             cached.value = some_expensive_function(T,P);
 *             cached.state1 = T;
 *             cached.state2 = P;
 *         }
 *         return cached.value;
 *     }
 * };
 * @endcode
 */
class ValueCache
{
public:
    //! Get a unique id for a cached value. Must be called exactly once for each
    //! method that implements caching behavior.
    int getId();

    //! Get a reference to a CachedValue object representing a scalar
    //! (double) with the given id.
    CachedScalar getScalar(int id) {
        return m_scalarCache[id];
    }

    //! Get a reference to a CachedValue object representing an array (vector<double>)
    //! with the given id.
    CachedArray getArray(int id) {
        return m_arrayCache[id];
    }

    //! Clear all cached values. This method should be called if the cached
    //! values may be invalidated in a way that is not represented by the state
    //! variables alone, such as a change to the constants defining a species
    //! thermodynamics as a function of temperature.
    void clear();

protected:
    //! Cached scalar values
    map<int, CachedValue<double>> m_scalarCache;

    //! Cached array values
    map<int, CachedValue<vector<double>>> m_arrayCache;

    //! The last assigned id. Automatically incremented by the getId() method.
    static int m_last_id;
};

}

#endif
