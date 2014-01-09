/*!
 * @file ct_thread.h
 * Header file containing utilities used to ensure thread safety.
 */

#ifndef CT_THREAD_H
#define CT_THREAD_H

#include "config.h"

#ifdef THREAD_SAFE_CANTERA
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>

#endif

namespace Cantera
{

#ifdef THREAD_SAFE_CANTERA

#if defined(BOOST_HAS_WINTHREADS)
typedef unsigned int cthreadId_t;
#elif defined(BOOST_HAS_PTHREADS)
typedef pthread_t cthreadId_t;
#endif

class thread_equal
{
public:
    bool operator()(cthreadId_t L, cthreadId_t R) {
#if defined(BOOST_HAS_WINTHREADS)
        return L == R;
#elif defined(BOOST_HAS_PTHREADS)
        return pthread_equal(L, R);
#endif
    }
};

typedef boost::mutex mutex_t;
typedef boost::mutex::scoped_lock ScopedLock;

#else
typedef int mutex_t;

class ScopedLock
{
public:
    explicit ScopedLock(const int m) : m_(m) {}
    int m_;
};

#endif // THREAD_SAFE_CANTERA

}

#endif
