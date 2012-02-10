/**
 * @file Storage.h
 */
#ifndef CTC_STORAGE_H
#define CTC_STORAGE_H

// Cantera includes
#include "kernel/Kinetics.h"
#include "kernel/TransportBase.h"

#include "Cabinet.h"
#include "clib_defs.h"


/**
 * Class to hold pointers to Cantera objects. Only one instance of
 * this class is needed.
 */
class Storage
{
public:
    Storage();
    virtual ~Storage();

    // vectors to hold pointers to objects
    std::vector<Cantera::Kinetics*> __ktable;
    std::vector<Cantera::thermo_t*> __thtable;
    std::vector<Cantera::Transport*> __trtable;

    std::map<std::string, int> __thmap;

    static Storage* storage() {
        if (__storage == 0) {
            __storage = new Storage;
        }
        return __storage;
    }


    size_t addThermo(Cantera::ThermoPhase* th);
    size_t addKinetics(Cantera::Kinetics* kin);
    size_t addTransport(Cantera::Transport* tr);
    //    int addNewTransport(int model, char* dbase, int th, int loglevel);
    int clear();
    void deleteKinetics(int n);
    void deleteThermo(int n);
    void deleteTransport(int n);
    size_t nThermo();
    static Storage* __storage;
};

inline Cantera::ThermoPhase* ph(int n)
{
    return Storage::__storage->__thtable[n];
}

inline Cantera::Kinetics* kin(int n)
{
    return Storage::__storage->__ktable[n];
}

inline Cantera::ThermoPhase* th(size_t n)
{
    return Storage::__storage->__thtable[n];
}

inline int thermo_index(std::string id)
{
    return Storage::__storage->__thmap[id];
}

inline Cantera::Transport* trans(int n)
{
    return Storage::__storage->__trtable[n];
}

#endif
