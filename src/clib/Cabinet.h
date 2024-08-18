/**
 * @file Cabinet.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CABINET_H
#define CT_CABINET_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/utilities.h"
#include <unordered_map>
#include <boost/range/adaptor/reversed.hpp>

namespace Cantera {

/**
 * Template for classes to hold pointers to objects. The SharedCabinet<M> class
 * maintains a list of pointers to objects of class M (or of subclasses of M). These
 * classes are used by the 'clib' interface library functions that provide access to
 * %Cantera C++ objects from outside C++. To refer to an existing object, the library
 * functions take an integer argument that specifies the location in the pointer list
 * maintained by the appropriate SharedCabinet<M> instance. The pointer is retrieved
 * from the list by the interface function, the desired method is invoked, and the
 * result returned to the non-C++ calling procedure. By storing the pointers in a
 * SharedCabinet, there is no need to encode them in a string or integer and pass
 * them out to the non-C++ calling routine, as some other interfacing schemes do.
 *
 * The SharedCabinet<M> class can be used to store pointers to arbitrary objects. In
 * most cases, class M is a base class with virtual methods, and the base class versions
 * of the methods throw CanteraError exceptions. The subclasses overload these methods
 * to implement the desired functionality. Class SharedCabinet<M> stores only the
 * base-class pointers, but since the methods are virtual, the method of the appropriate
 * subclass will be invoked.
 *
 * As the SharedCabinet<M> class uses smart pointers, it is set up to allow deleting
 * objects in an inherently safe manner. Method 'del' does the following. If called
 * with n >= 0, it dereferences the object. The original object is only destroyed if the
 * reference is not shared by other objects. In this way, if it is deleted again
 * inadvertently nothing happens, and if an attempt is made to reference the object by
 * its index number, a standard exception is thrown.
 *
 * The SharedCabinet<M> class is implemented as a singleton. The constructor is never
 * explicitly called; instead, static function SharedCabinet<M>::SharedCabinet() is
 * called to obtain a pointer to the instance. This function calls the constructor on
 * the first call and stores the pointer to this instance. Subsequent calls simply
 * return the already-created pointer.
 */
template<class M>
class SharedCabinet
{
public:
    typedef vector<shared_ptr<M>>& dataRef;
    typedef std::unordered_map<const M*, set<int>>& lookupRef;

    /**
     * Constructor.
     */
    SharedCabinet() {}

    /**
     * Add a new object. The index of the object is returned.
     */
    static int add(shared_ptr<M> obj, int parent=-1) {
        dataRef data = getData();
        data.push_back(obj);
        auto& parents = getParents();
        parents.push_back(parent);
        int idx = static_cast<int>(data.size()) - 1;
        lookupRef lookup = getLookup();
        if (lookup.count(obj.get())) {
            lookup[obj.get()].insert(idx);
        } else {
            lookup[obj.get()] = {idx};
        }
        return idx;
    }

    /**
     * Return cabinet size.
     */
    static int size() {
        return static_cast<int>(getData().size());
    }

    /**
     * Delete all objects without erasing mapping.
     */
    static int clear() {
        dataRef data = getData();
        for (size_t i = 0; i < data.size(); i++) {
            del(static_cast<int>(i));
        }
        return 0;
    }

    /**
     * Delete all objects and erase mapping.
     */
    static int reset() {
        getData().clear();
        getParents().clear();
        getLookup().clear();
        return 0;
    }

    /**
     * Add a copy of the nth object to storage. The index of the new entry is returned.
     */
    static int copy(int n) {
        dataRef data = getData();
        try {
            return add(*data[n]);  // do not copy parent to avoid ambiguous data
        } catch (std::exception& err) {
            throw CanteraError("SharedCabinet::newCopy", err.what());
        }
    }

    /**
     * Delete the nth object.
     */
    static void del(int n) {
        dataRef data = getData();
        if (n >= 0 && n < len(data)) {
            lookupRef lookup = getLookup();
            if (!lookup.count(data[n].get())) {
                throw CanteraError("SharedCabinet::del",
                    "Lookup table does not contain reference to object.");
            }
            if (lookup[data[n].get()].size() == 1) {
                // set only contains one index
                lookup.erase(data[n].get());
            } else {
                // remove index n from the reverse lookup table
                lookup[data[n].get()].erase(n);
            }
            data[n].reset();
        } else {
            throw CanteraError("SharedCabinet::del",
                "Attempt made to delete a non-existing object.");
        }
    }

    /**
     * Return handle of parent to object n.
     */
    static int parent(int n) {
        auto& parents = getParents();
        if (n < 0 || n >= len(parents)) {
            throw CanteraError("SharedCabinet::parent", "Index {} out of range.", n);
        }
        return parents[n];
    }

    /**
     * Return a shared pointer to object n.
     */
    static shared_ptr<M>& at(int n) {
        dataRef data = getData();
        if (n < 0 || n >= len(data)) {
            throw CanteraError("SharedCabinet::at", "Index {} out of range.", n);
        }
        if (!data[n]) {
            throw CanteraError("SharedCabinet::at",
                "Object with index {} has been deleted.", n);
        }
        return data[n];
    }

    /**
     * Return object n, cast to the specified type.
     */
    template <class T>
    static shared_ptr<T> as(int n) {
        auto obj = std::dynamic_pointer_cast<T>(at(n));
        if (obj) {
            return obj;
        }
        throw CanteraError("SharedCabinet::as", "Item is not of the correct type.");
    }

    /**
     * Return a reference to object n, cast to a reference of the specified type.
     */
    template <class T>
    static T& get(int n) {
        auto obj = std::dynamic_pointer_cast<T>(at(n));
        if (obj) {
            return *obj;
        }
        throw CanteraError("SharedCabinet::get", "Item is not of the correct type.");
    }

    /**
     * Return the index in the SharedCabinet to the specified object, or -1 if the
     * object is not in the SharedCabinet. If multiple indices reference the same
     * object, the index of the last one added is returned.
     */
    static int index(const M& obj, int parent=-1) {
        lookupRef lookup = getLookup();
        if (!lookup.count(&obj)) {
            return -1;
        }
        set<int>& entry = lookup.at(&obj);
        auto& parents = getParents();
        for (const auto e : boost::adaptors::reverse(entry)) {
            if (parents[e] == parent) {
                return e;
            }
        }
        return -2;  // not found
    }

private:
    /**
     * Static function that returns a pointer to the data member of
     * the singleton SharedCabinet<M> instance. All member functions should
     * access the data through this function.
     */
    static dataRef getData() {
        if (s_storage == nullptr) {
            s_storage = new SharedCabinet<M>();
        }
        return s_storage->m_table;
    }

    /**
     * Static function that returns a pointer to the list of parent object handles of
     * the singleton SharedCabinet<M> instance. All member functions should
     * access the data through this function.
     */
    static vector<int>& getParents() {
        if (s_storage == nullptr) {
            s_storage = new SharedCabinet<M>();
        }
        return s_storage->m_parents;
    }

    /**
     * Static function that returns a pointer to the reverse lookup table of
     * the singleton SharedCabinet<M> instance. All member functions should
     * access the lookup table through this function.
     */
    static lookupRef getLookup() {
        if (s_storage == nullptr) {
            s_storage = new SharedCabinet<M>();
        }
        return s_storage->m_lookup;
    }

    /**
     * Pointer to the single instance of this class.
     */
    static SharedCabinet<M>* s_storage;

    /**
     * Reverse lookup table for the single instance of this class.
     */
    std::unordered_map<const M*, set<int>> m_lookup;

    /**
     * List to hold handles of parent objects.
     */
    vector<int> m_parents;

    /**
     * List to hold pointers to objects.
     */
    vector<shared_ptr<M>> m_table;
};

}

#endif
