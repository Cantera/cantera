//! @file global.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/FactoryBase.h"
#include "application.h"
#include "cantera/base/AnyMap.h"

#define BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED
#include <boost/stacktrace.hpp>
#include <boost/core/demangle.hpp>

#include <signal.h>

namespace Cantera
{

//! Return a pointer to the application object
static Application* app()
{
    return Application::Instance();
}

// **************** Text Logging ****************

void setLogger(Logger* logwriter)
{
    try {
        app()->setLogger(logwriter);
    } catch (const std::bad_alloc&) {
        logwriter->error("bad alloc thrown by app()");
    }
}

void setLogger(unique_ptr<Logger> logwriter)
{
    try {
        app()->setLogger(std::move(logwriter));
    } catch (const std::bad_alloc&) {
        logwriter->error("bad alloc thrown by app()");
    }
}

void writelog_direct(const string& msg)
{
    app()->writelog(msg);
}

void writelogendl()
{
    app()->writelogendl();
}

void writeline(char repeat, size_t count, bool endl_after, bool endl_before)
{
    if (endl_before) {
        writelogendl();
    }
    writelog_direct(string(count, repeat));
    if (endl_after) {
        writelogendl();
    }
}

void _warn_deprecated(const string& method, const string& extra)
{
    app()->warn_deprecated(method, extra);
}

void _warn(const string& warning, const string& method, const string& extra)
{
    app()->warn(warning, method, extra);
}

void suppress_deprecation_warnings()
{
    app()->suppress_deprecation_warnings();
}

void make_deprecation_warnings_fatal()
{
    app()->make_deprecation_warnings_fatal();
}

void suppress_warnings()
{
    app()->suppress_warnings();
}

bool warnings_suppressed()
{
    return app()->warnings_suppressed();
}

void make_warnings_fatal()
{
    app()->make_warnings_fatal();
}

void suppress_thermo_warnings(bool suppress)
{
    app()->suppress_thermo_warnings(suppress);
}

bool thermo_warnings_suppressed()
{
    return app()->thermo_warnings_suppressed();
}

void use_legacy_rate_constants(bool legacy)
{
    app()->use_legacy_rate_constants(legacy);
}

bool legacy_rate_constants_used()
{
    return app()->legacy_rate_constants_used();
}

namespace {
void stacktraceWriter(int signum) {
    if (signum == SIGSEGV) {
        std::cerr << "\nSegmentation fault. Stack trace:\n";
    } else {
        std::cerr << "\nProcess terminated abnormally. Stack trace:\n";
    }
    ::signal(signum, SIG_DFL);
    std::cerr << boost::stacktrace::stacktrace();
    ::raise(signum);
}
}

void printStackTraceOnSegfault()
{
    // Install signal handler to print stacktrace in case of segfault.
    ::signal(SIGSEGV, &stacktraceWriter);
    ::signal(SIGABRT, &stacktraceWriter);
}

// **************** Global Data ****************

void appdelete()
{
    Application::ApplicationDestroy();
    FactoryBase::deleteFactories();
}

void thread_complete()
{
    app()->thread_complete();
}

string version()
{
    return CANTERA_VERSION;
}

string gitCommit()
{
#ifdef GIT_COMMIT
    return GIT_COMMIT;
#else
    return "unknown";
#endif
}

void addDirectory(const string& dir)
{
    app()->addDataDirectory(dir);
}

string getDataDirectories(const string& sep)
{
    return app()->getDataDirectories(sep);
}

string findInputFile(const string& name)
{
    return app()->findInputFile(name);
}

void loadExtension(const string& extType, const string& name)
{
    app()->loadExtension(extType, name);
}

void loadExtensions(const AnyMap& node)
{
    if (!node.hasKey("extensions")) {
        return;
    }
    for (auto& extension : node["extensions"].asVector<AnyMap>()) {
        loadExtension(extension["type"].asString(), extension["name"].asString());
    }
}

void searchPythonVersions(const string& versions) {
    app()->searchPythonVersions(versions);
}

bool debugModeEnabled()
{
#ifdef NDEBUG
    return false;
#else
    return true;
#endif
}

bool usesHDF5()
{
#if CT_USE_HDF5
    return true;
#else
    return false;
#endif
}

vector<FactoryBase*> FactoryBase::s_vFactoryRegistry;

string demangle(const std::type_info& type)
{
    static map<string, string> typenames = {
        {typeid(void).name(), "void"},
        {typeid(double).name(), "double"},
        {typeid(long int).name(), "long int"},
        {typeid(bool).name(), "bool"},
        {typeid(string).name(), "string"},
        {typeid(vector<AnyValue>).name(), "vector<AnyValue>"},
        {typeid(vector<AnyMap>).name(), "vector<AnyMap>"},
        {typeid(vector<double>).name(), "vector<double>"},
        {typeid(vector<long int>).name(), "vector<long int>"},
        {typeid(vector<bool>).name(), "vector<bool>"},
        {typeid(vector<string>).name(), "vector<string>"},
        {typeid(vector<vector<double>>).name(), "vector<vector<double>>"},
        {typeid(vector<vector<long int>>).name(), "vector<vector<long int>>"},
        {typeid(vector<vector<bool>>).name(), "vector<vector<bool>>"},
        {typeid(vector<vector<string>>).name(), "vector<vector<string>>"},
        {typeid(AnyMap).name(), "AnyMap"},
    };

    if (typenames.count(type.name())) {
        return typenames[type.name()];
    } else {
        return boost::core::demangle(type.name());
    }
}

}
