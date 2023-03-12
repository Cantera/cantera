/**
 * @file Cabinet.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CABINET_H
#define CT_CABINET_H

#include "cantera/base/ctexceptions.h"
#include <unordered_map>

namespace Cantera {

/**
 * Template for classes to hold pointers to objects. The Cabinet<M>
 * class maintains a list of pointers to objects of class M (or of
 * subclasses of M).  These classes are used by the 'clib' interface
 * library functions that provide access to Cantera C++ objects from
 * outside C++. To refer to an existing object, the library functions
 * take an integer argument that specifies the location in the pointer
 * list maintained by the appropriate Cabinet<M> instance. The pointer
 * is retrieved from the list by the interface function, the desired
 * method is invoked, and the result returned to the non-C++ calling
 * procedure. By storing the pointers in a 'cabinet', there is no need
 * to encode them in a std::string or integer and pass them out to the
 * non-C++ calling routine, as some other interfacing schemes do.
 *
 * The Cabinet<M> class can be used to store pointers to any class
 * that is default-constructible (that is, has a constructor that takes
 * no arguments). The requirement that the class be
 * default-constructible arises since the Cabinet constructor always
 * creates an instance of M by invoking 'new M', and stores a pointer
 * to it as the first entry in the list. In most cases, class M is a
 * base class with virtual methods, and the base class versions of the
 * methods throw CanteraError exceptions. The subclasses overload
 * these methods to implement the desired functionality. Class
 * Cabinet<M> stores only the base-class pointers, but since the
 * methods are virtual, the method of the appropriate subclass will be
 * invoked.
 *
 * The Cabinet<M> class is set up to allow deleting objects in a safe
 * manner, *provided* that method 'delete' is used, and the destructor
 * for the object is not called directly. Method 'delete' does the
 * following. If called with n = 0, it does nothing, since the first
 * object in the list (the default-constructed base class instance) is
 * never destroyed. If called with n > 0, it deletes the object, and
 * replaces the pointer to where the object had been (but is no more)
 * with a pointer to the first object. In this way, if it is deleted
 * again inadvertently nothing happens, and if an attempt is made to
 * reference the object by its index number, the base-class object
 * will be referenced instead, which will throw an exception.  If
 * instead the pointer were stored in the referring code, there would
 * always be the chance that
 *
 * The Cabinet<M> class is implemented as a singlet. The constructor
 * is never explicitly called; instead, static function
 * Cabinet<M>::cabinet() is called to obtain a pointer to the
 * instance. This function calls the constructor on the first call and
 * stores the pointer to this instance. Subsequent calls simply return
 * the already-created pointer.
 *
 * Set canDelete to false if the 'clear' method should not delete the entries.
 */

template<class M, bool canDelete=true>
class Cabinet
{
public:
    typedef std::vector<M*>& dataRef;
    /**
     * Destructor. Delete all objects in the list.
     */
    virtual ~Cabinet() {
        clear();
    }

    /**
     * Add a new object. The index of the object is returned.
     */
    static int add(M* ptr) {
        dataRef data = getData();
        data.push_back(ptr);
        return static_cast<int>(data.size()) - 1;
    }

    /**
     * Make a new copy of an existing object.  The index of the new
     * object is returned.
     */
    static int newCopy(int i) {
        dataRef data = getData();
        try {
            M* old = data[i];
            data.push_back(new M(*old));
            return static_cast<int>(data.size()) - 1;
        } catch (std::exception& err) {
            throw Cantera::CanteraError("Cabinet::newCopy", err.what());
        }
    }

    /**
     * Delete all objects but the first.
     */
    static int clear() {
        dataRef data = getData();
        for (size_t i = 1; i < data.size(); i++) {
            del(i);
        }
        if (canDelete) {
            delete data[0];
        }
        data.clear();
        add(new M);
        return 0;
    }

    /**
     * Delete the nth object. After the object is deleted, the pointer
     * to it in the list is replaced by a pointer to the first element
     * in the list.
     */
    static void del(size_t n) {
        dataRef data = getData();
        if (n == 0) {
            return;
        }
        if (data[n] != data[0]) {
            if (canDelete) {
                delete data[n];
            }
            data[n] = data[0];
        } else {
            throw Cantera::CanteraError("Cabinet::del",
                                        "Attempt made to delete an already-deleted object.");
        }
    }

    /**
     * Return a reference to object n.
     */
    static M& item(size_t n) {
        dataRef data = getData();
        if (n < data.size()) {
            return *data[n];
        } else {
            throw Cantera::CanteraError("Cabinet::item","index out of range {}", n);
        }
    }

    /**
     * Return a reference to object n, cast to a reference of the specified type.
     */
    template <class T>
    static T& get(size_t n) {
        T* x = dynamic_cast<T*>(&item(n));
        if (x == 0) {
            throw Cantera::CanteraError("Cabinet::get",
                                        "Item is not of the correct type.");
        }
        return *x;
    }

    /**
     * Return the index in the Cabinet to the specified object, or -1
     * if the object is not in the cabinet.
     */
    static int index(const M& obj) {
        dataRef data = getData();
        auto loc = std::find(data.begin(), data.end(), &obj);
        if (loc != data.end()) {
            return static_cast<int>(loc-data.begin());
        } else {
            return -1;
        }
    }

    /**
     * Constructor.
     */
    Cabinet() {
        m_table.push_back(new M);
    }

private:
    /**
     * Static function that returns a pointer to the data member of
     * the singleton Cabinet<M> instance. All member functions should
     * access the data through this function.
     */
    static dataRef getData() {
        if (s_storage == 0) {
            s_storage = new Cabinet<M, canDelete>();
        }
        return s_storage->m_table;
    }

    /**
     * Pointer to the single instance of this class.
     */
    static Cabinet<M, canDelete>* s_storage;

    /**
     * Vector to hold pointers to objects.
     */
    std::vector<M*> m_table;
};


/**
 * Template for classes to hold pointers to objects. The SharedCabinet<M> class
 * maintains a list of pointers to objects of class M (or of subclasses of M). These
 * classes are used by the 'clib' interface library functions that provide access to
 * Cantera C++ objects from outside C++. To refer to an existing object, the library
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
 * The SharedCabinet<M> class is implemented as a singlet. The constructor is never
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
    static int add(shared_ptr<M> obj) {
        dataRef data = getData();
        data.push_back(obj);
        int idx = data.size() - 1;
        lookupRef lookup = getLookup();
        if (index(*obj) >= 0) {
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
        int size = getData().size();
        return size;
    }

    /**
     * Delete all objects without erasing mapping.
     */
    static int clear() {
        dataRef data = getData();
        for (size_t i = 0; i < data.size(); i++) {
            del(i);
        }
        return 0;
    }

    /**
     * Delete all objects and erase mapping.
     */
    static int reset() {
        getData().clear();
        getLookup().clear();
        return 0;
    }

    /**
     * Add a copy of the nth object to storage. The index of the new entry is returned.
     */
    static int copy(size_t n) {
        dataRef data = getData();
        try {
            return add(*data[n]);
        } catch (std::exception& err) {
            throw CanteraError("SharedCabinet::newCopy", err.what());
        }
    }

    /**
     * Delete the nth object.
     */
    static void del(size_t n) {
        dataRef data = getData();
        if (n >= 0 && n < data.size()) {
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
     * Return a shared pointer to object n.
     */
    static shared_ptr<M>& at(size_t n) {
        dataRef data = getData();
        if (n < 0 || n >= data.size()) {
            throw CanteraError("SharedCabinet::at", "Index {} out of range.", n);
        }
        return data[n];
    }

    /**
     * Return object n, cast to the specified type.
     */
    template <class T>
    static shared_ptr<T> as(size_t n) {
        auto obj = std::dynamic_pointer_cast<T>(at(n));
        if (obj) {
            return obj;
        }
        throw CanteraError("SharedCabinet::as", "Item is not of the correct type.");
    }

    /**
     * Return a reference to object n.
     */
    static M& item(size_t n) {
        auto ptr = at(n);
        if (!ptr) {
            throw CanteraError("SharedCabinet::item",
                "Object with index {} has been deleted.", n);
        }
        return *ptr;
    }

    /**
     * Return a reference to object n, cast to a reference of the specified type.
     */
    template <class T>
    static T& get(size_t n) {
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
    static int index(const M& obj) {
        lookupRef lookup = getLookup();
        if (!lookup.count(&obj)) {
            return -1;
        }
        set<int>& entry = lookup.at(&obj);
        return *entry.rbegin();
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
     * list to hold pointers to objects.
     */
    vector<shared_ptr<M>> m_table;
};

}

#endif
