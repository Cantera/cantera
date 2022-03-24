//! @file global.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/FactoryBase.h"
#include "cantera/base/xml.h"
#include "application.h"
#include "units.h"
#include "cantera/base/AnyMap.h"
#ifdef CT_USE_DEMANGLE
  #include <boost/core/demangle.hpp>
#endif

using namespace std;

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

void writelog_direct(const std::string& msg)
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
    writelog_direct(std::string(count, repeat));
    if (endl_after) {
        writelogendl();
    }
}

void _warn_deprecated(const std::string& method, const std::string& extra)
{
    app()->warn_deprecated(method, extra);
}

void _warn(const std::string& warning,
           const std::string& method, const std::string& extra)
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

// **************** Global Data ****************

Unit* Unit::s_u = 0;
std::mutex Unit::units_mutex;

void appdelete()
{
    Application::ApplicationDestroy();
    FactoryBase::deleteFactories();
    Unit::deleteUnit();
}

void thread_complete()
{
    app()->thread_complete();
}

std::string gitCommit()
{
#ifdef GIT_COMMIT
    return GIT_COMMIT;
#else
    return "unknown";
#endif
}

XML_Node* get_XML_File(const std::string& file, int debug)
{
    return app()->get_XML_File(file, debug);
}

XML_Node* get_XML_from_string(const std::string& text)
{
    return app()->get_XML_from_string(text);
}

void close_XML_File(const std::string& file)
{
    app()->close_XML_File(file);
}

void addDirectory(const std::string& dir)
{
    app()->addDataDirectory(dir);
}

std::string getDataDirectories(const std::string& sep)
{
    return app()->getDataDirectories(sep);
}

std::string findInputFile(const std::string& name)
{
    return app()->findInputFile(name);
}

doublereal toSI(const std::string& unit)
{
    doublereal f = Unit::units()->toSI(unit);
    if (f) {
        return f;
    } else {
        throw CanteraError("toSI","unknown unit string: "+unit);
    }
    return 1.0;
}

doublereal actEnergyToSI(const std::string& unit)
{
    doublereal f = Unit::units()->actEnergyToSI(unit);
    if (f) {
        return f;
    }
    return 1.0;
}

string canteraRoot()
{
    warn_deprecated("canteraRoot",
                    "Unused in Cantera. To be removed after Cantera 2.6");
    char* ctroot = getenv("CANTERA_ROOT");
    if (ctroot != 0) {
        return string(ctroot);
    }
#ifdef CANTERA_ROOT
    return string(CANTERA_ROOT);
#else
    return "";
#endif

}

bool debugModeEnabled()
{
#ifdef NDEBUG
    return false;
#else
    return true;
#endif
}

//! split a string at a '#' sign. Used to separate a file name from an id string.
/*!
 *   @param    src     Original string to be split up. This is unchanged.
 *   @param    file    Output string representing the first part of the string,
 *       which is the filename.
 *   @param    id      Output string representing the last part of the string,
 *       which is the id.
 */
static void split_at_pound(const std::string& src, std::string& file, std::string& id)
{
    string::size_type ipound = src.find('#');
    if (ipound != string::npos) {
        id = src.substr(ipound+1,src.size());
        file = src.substr(0,ipound);
    } else {
        id = "";
        file = src;
    }
}

XML_Node* get_XML_Node(const std::string& file_ID, XML_Node* root)
{
    std::string fname, idstr;
    XML_Node* db;
    split_at_pound(file_ID, fname, idstr);
    if (fname == "") {
        if (!root) throw CanteraError("get_XML_Node",
                                          "no file name given. file_ID = "+file_ID);
        db = root->findID(idstr, 3);
    } else {
        try {
            findInputFile(fname);
        } catch (CanteraError& err) {
            // See if the input file can be found with a different format
            if (fname.rfind(".xml") == fname.size() - 4) {
                fname.replace(fname.size() - 3, 3, "cti");
            } else if (fname.rfind(".cti") == fname.size() - 4) {
                fname.replace(fname.size() - 3, 3, "xml");
            }
            try {
                findInputFile(fname);
            } catch (CanteraError&) {
                // rethrow the original error, which indicates the given file name
                throw err;
            }
        }
        XML_Node* doc = get_XML_File(fname);
        if (!doc) throw CanteraError("get_XML_Node",
                                         "get_XML_File failed trying to open "+fname);
        db = doc->findID(idstr, 3);
    }
    if (!db) {
        throw CanteraError("get_XML_Node",
                           "id tag '"+idstr+"' not found.");
    }
    return db;
}

XML_Node* get_XML_NameID(const std::string& nameTarget,
                         const std::string& file_ID,
                         XML_Node* root)
{
    string fname, idTarget;
    XML_Node* db;
    split_at_pound(file_ID, fname, idTarget);
    if (fname == "") {
        if (!root) {
            return 0;
        }
        db = root->findNameID(nameTarget, idTarget);
    } else {
        XML_Node* doc = get_XML_File(fname);
        if (!doc) {
            return 0;
        }
        db = doc->findNameID(nameTarget, idTarget);
    }
    return db;
}

std::vector<FactoryBase*> FactoryBase::s_vFactoryRegistry;

std::string demangle(const std::type_info& type)
{
    static std::map<std::string, std::string> typenames = {
        {typeid(void).name(), "void"},
        {typeid(double).name(), "double"},
        {typeid(long int).name(), "long int"},
        {typeid(bool).name(), "bool"},
        {typeid(std::string).name(), "string"},
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
        #ifdef CT_USE_DEMANGLE
            return boost::core::demangle(type.name());
        #else
            return type.name();
        #endif
    }
}

}
