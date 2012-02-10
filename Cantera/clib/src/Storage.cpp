/**
 * @file Storage.cpp
 */
// Cantera includes
#include "kernel/Kinetics.h"
#include "kernel/TransportFactory.h"

#include "Storage.h"

using namespace std;
using namespace Cantera;


Storage::Storage()
{
    addThermo(new ThermoPhase);
    addKinetics(new Kinetics);
    addTransport(newTransportMgr());
}

Storage::~Storage()
{
    clear();
}

size_t Storage::addThermo(thermo_t* th)
{
    if (th->index() != npos) {
        return th->index();
    }
    __thtable.push_back(th);
    size_t n = __thtable.size() - 1;
    th->setIndex(n);
    //string id = th->id();
    //if (__thmap.count(id) == 0) {
    //    __thmap[id] = n;
    //    th->setID(id);
    //}
    //else {
    //    throw CanteraError("Storage::addThermo","id already used");
    //    return -1;
    //}
    return n;
}

size_t Storage::nThermo()
{
    return __thtable.size();
}

size_t Storage::addKinetics(Kinetics* kin)
{
    if (kin->index() != npos) {
        return kin->index();
    }
    __ktable.push_back(kin);
    size_t n = __ktable.size() - 1;
    kin->setIndex(n);
    return n;
}

size_t Storage::addTransport(Transport* tr)
{
    if (tr->index() != npos) {
        return tr->index();
    }
    __trtable.push_back(tr);
    size_t n = __trtable.size() - 1;
    tr->setIndex(n);
    return n;
}

// int Storage::addNewTransport(int model, char* dbase, int th,
//     int loglevel) {
//     try {
//         ThermoPhase* thrm = __thtable[th];
//         Transport* tr = newTransportMgr(model,
//             string(dbase), thrm, loglevel);
//         __trtable.push_back(tr);
//         return __trtable.size() - 1;
//     }
//     catch (CanteraError) {return -1;}
//     catch (...) {return ERR;}
// }

int Storage::clear()
{
    size_t i, n;
    n = __thtable.size();
    for (i = 1; i < n; i++) {
        if (__thtable[i] != __thtable[0]) {
            delete __thtable[i];
            __thtable[i] = __thtable[0];
        }
    }
    n = __ktable.size();
    for (i = 1; i < n; i++) {
        if (__ktable[i] != __ktable[0]) {
            delete __ktable[i];
            __ktable[i] = __ktable[0];
        }
    }
    n = __trtable.size();
    for (i = 1; i < n; i++) {
        if (__trtable[i] != __trtable[0]) {
            delete __trtable[i];
            __trtable[i] = __trtable[0];
        }
    }
    return 0;
}

void Storage::deleteKinetics(int n)
{
    if (n == 0) {
        return;
    }
    if (__ktable[n] != __ktable[0]) {
        delete __ktable[n];
    }
    __ktable[n] = __ktable[0];
}

void Storage::deleteThermo(int n)
{
    if (n == 0) {
        return;
    }
    if (n < 0 || n >= (int) __thtable.size()) {
        throw CanteraError("deleteThermo","illegal index");
    }
    __thtable[n] = __thtable[0];
}

void Storage::deleteTransport(int n)
{
    if (n == 0) {
        return;
    }
    if (__trtable[n] != __trtable[0]) {
        delete __trtable[n];
    }
    __trtable[n] = __trtable[0];
}

Storage* Storage::__storage = 0;

