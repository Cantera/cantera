#include "cantera/base/ctexceptions.h"

#include "application.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

namespace Cantera {

// *** Exceptions ***

CanteraError::CanteraError(std::string procedure, std::string msg)
{
    Application::Instance()->addError(procedure, msg);
}

ArraySizeError::ArraySizeError(std::string proc, size_t sz, size_t reqd) :
    CanteraError(proc, "Array size ("+int2str(sz)+
                 ") too small. Must be at least "+int2str(reqd))
{
}


ElementRangeError::ElementRangeError(std::string func, size_t m, size_t mmax) :
    CanteraError(func, "Element index " + int2str(m) +
                 " outside valid range of 0 to " + int2str(mmax-1))
{
}

// *** Warnings ***

void deprecatedMethod(std::string classnm, std::string oldnm, std::string newnm)
{
    writelog(">>>> WARNING: method "+oldnm+" of class "+classnm
             +" is deprecated.\n");
    writelog("         Use method "+newnm+" instead.\n");
    writelog("         (If you want to rescue this method from deprecated\n");
    writelog("         status, see http://www.cantera.org/deprecated.html)");
}

void removeAtVersion(std::string func, std::string version)
{
    writelog("Removed procedure: "+func+"\n");
    writelog("Removed in version: "+version+"\n");
    throw CanteraError("removeAtVersion: "+ func,"procedure has been removed.");
}


} // namespace Cantera
