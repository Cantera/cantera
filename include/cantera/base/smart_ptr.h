#ifndef CT_SMART_PTR
#define CT_SMART_PTR

#include "config.h"

#if defined CT_USE_STD_SHARED_PTR
#include <memory>
namespace Cantera
{
    using std::shared_ptr;
}

#elif defined CT_USE_TR1_SHARED_PTR
#include <tr1/memory>
namespace Cantera
{
    using std::tr1::shared_ptr;
}

#elif defined CT_USE_MSFT_SHARED_PTR
#include <memory>
namespace Cantera
{
    using std::tr1::shared_ptr;
}

#elif defined CT_USE_BOOST_SHARED_PTR
#include <boost/shared_ptr.hpp>
namespace Cantera
{
    using boost::shared_ptr;
}

#else
#error "No shared_ptr implementation available"
#endif

#endif
