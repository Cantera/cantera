/**
 * @file Cabinet.h
 */
/*
 *      $Id: Cabinet.h,v 1.9 2009/07/11 17:16:09 hkmoffa Exp $
 */


#ifndef CT_CABINET_H
#define CT_CABINET_H

#include <vector>
#include "stringUtils.h"
#include "config.h"

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
 * instead the pointer were stored in the refering code, there would
 * always be the chance that
 *
 * The Cabinet<M> class is implemented as a singlet. The constructor
 * is never explicitly called; instead, static function
 * Cabinet<M>::cabinet() is called to obtain a pointer to the
 * instance. This function calls the constructor on the first call and
 * stores the pointer to this instance. Subsequent calls simply return
 * the already-created pointer.
 */

template<class M>
class Cabinet {
public:

    /**
     * Destructor. Delete all objects in the list.
     */
    virtual ~Cabinet() {clear();}

    /**
     * Static function that returns a pointer to the one Cabinet<M>
     * instance. All access to the Cabinet<M> instance should go
     * through this function.
     */
    static Cabinet<M>* cabinet(bool canDelete = true) {
        if (__storage == 0) {
            __storage = new Cabinet<M>(canDelete);
        }
        return __storage;
    }


    /** 
     * Add a new object. The index of the object is returned.  
     */
    int add(M* ptr) {
        //try {
        __table.push_back(ptr);
        return static_cast<int>(__table.size()) - 1;
        //}
        //catch (CanteraError) {return -1;}
        //catch (...) {return -999;}
    }


    /** 
     * Make a new copy of an existing object.  The index of the new
     * object is returned.
     */
    int newCopy(int i) {
        try {
            M* old = __table[i];
            __table.push_back(new M(*old));
            return static_cast<int>(__table.size()) - 1;
        }
        catch (Cantera::CanteraError) {return -1;}
        catch (...) {return -999;}
    }


    /** 
     * Assign one object (index j) to another (index i).  This method
     * is not used currently, and may be removed from the class in the
     * future.
     */
    int assign(int i, int j) {
        try {
            M* src = __table[j];
            M* dest = __table[i];
            *dest = *src;
            return 0;
        }
        catch (Cantera::CanteraError) {return -1;}
        catch (...) {return -999;}
    }


    /** 
     * Delete all objects but the first.
     */
    int clear() {
        int i, n;
        n = static_cast<int>(__table.size());
        for (i = 1; i < n; i++) {del(i);}
        if (_can_delete) delete __table[0];
        __table.clear();
        add(new M);
        return 0;
    }


    /**
     * Delete the nth object. After the object is deleted, the pointer
     * to it in the list is replaced by a pointer to the first element
     * in the list.
     */
    void del(int n) {
        if (n == 0) return;
        if (__table[n] != __table[0]) {
            if (_can_delete) delete __table[n];
            __table[n] = __table[0]; 
        }
        else {
            throw Cantera::CanteraError("Cabinet<M>::del", 
                "Attempt made to delete an already-deleted object.");
        } 
    }


    /** 
     * Return a pointer to object n.
     */
    M* item(int n) {
        if (n >= 0 && n < int(__table.size()))
            return __table[n];
        else {
            throw Cantera::CanteraError("item","index out of range"+Cantera::int2str(n));
            //return __table[0];
        }
    }

    /**
     * Constructor. 
     */
    Cabinet(bool canDelete = true) : _can_delete(canDelete) { add(new M); }

private:

    /**
     * Constructor. 
     */
    //    Cabinet(bool canDelete = true) : _can_delete(canDelete) { add(new M); }


    /**
     * Pointer to the single instance of this class.
     */
    static Cabinet<M>* __storage;

    /**
     * Vector to hold pointers to objects.
     */
    std::vector<M*> __table;

    /**
     * Set to false if 'clear' should not delete the entries.
     */
    bool _can_delete;
};

//! Declaration stating that the storage for the static member
//! of each instanteated template will exist
/*!
 *   The actual storage will be allocated in .cpp files
 */
#ifdef NEEDS_GENERIC_TEMPL_STATIC_DECL
template<class M>  Cabinet<M>*       Cabinet<M>::__storage;
#endif

#endif
