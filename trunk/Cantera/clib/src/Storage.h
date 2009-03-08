#ifndef CTC_STORAGE_H
#define CTC_STORAGE_H

// Cantera includes
#include "Kinetics.h"
#include "transport/TransportBase.h"

#include "Cabinet.h"
#include "clib_defs.h"


/**
 * Class to hold pointers to Cantera objects. Only one instance of
 * this class is needed.
 */
class Storage {
public:
    Storage();
    virtual ~Storage();

    // vectors to hold pointers to objects
    vector<Kinetics*> __ktable;
    vector<thermo_t*> __thtable;
    vector<Transport*> __trtable;

    map<string, int> __thmap;

    static Storage* storage() {
        if (__storage == 0) {
            __storage = new Storage;
        }
        return __storage;
    }


    int addThermo(ThermoPhase* th);
    int addKinetics(Kinetics* kin);
    int addTransport(Transport* tr);
    //    int addNewTransport(int model, char* dbase, int th, int loglevel);
    int clear();
    void deleteKinetics(int n);
    void deleteThermo(int n);
    void deleteTransport(int n);
    int nThermo();
    static Storage* __storage;
};

inline ThermoPhase* ph(int n) {
    return Storage::__storage->__thtable[n];
}

inline Kinetics* kin(int n) {
    return Storage::__storage->__ktable[n];
}

inline ThermoPhase* th(int n) {
    return Storage::__storage->__thtable[n];
}

inline int thermo_index(string id) {
    return Storage::__storage->__thmap[id];
}

inline Transport* trans(int n) {
    return Storage::__storage->__trtable[n];
}

#endif
