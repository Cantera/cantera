/**
 * @file Cabinet.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CABINET_H
#define CT_CABINET_H

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"
#include "clib_utils.h"

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
 * that is default-constructible (i.e., has a constructor that takes
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
        } catch (...) {
            return Cantera::handleAllExceptions(-1, -999);
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

//! Declaration stating that the storage for the static member
//! of each instantiated template will exist
/*!
 *   The actual storage will be allocated in .cpp files
 */
#ifdef NEEDS_GENERIC_TEMPL_STATIC_DECL
template<class M, bool canDelete> Cabinet<M, canDelete>* Cabinet<M, canDelete>::s_storage;
#endif

#endif
