/**
 *  @file PropertyUpdater.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

TO BE REMOVED


// Copyright 2001  California Institute of Technology


#ifndef CT_PROPUPDATER_H
#define CT_PROPUPDATER_H

#include "ct_defs.h"

/*
 * Reacting flow simulations require the evaluation of many quantities
 * that are expensive to compute, such as reaction rates of progress,
 * equilibrium constants, and multicomponent transport properties. In
 * many cases, these quantities in turn require evaluating properties
 * that may depend only on temperature, composition, or other
 * something else. For example, reaction rate coefficients depend only
 * on temperature (or may be decomposed into parts that
 * do). Re-evaluating the reaction rate coefficients each time the
 * rates of progress of the reactions are computed is inefficient,
 * since the temperature may not have changed since the last
 * call. This occurs commonly when using implicit methods that require
 * evaluating a Jacobian, and the Jacobian is evaluated numerically.
 *
 * Cantera implements a technique to manage property updating that
 * results in properties being re-evaluated only when the quantities
 * they depend on have changed, and only when the properties are
 * needed. 
 *
 * The basic idea is as follows. For every property that needs
 * updating, a class is derived from class Updater, and its method
 * 'update()' is overloaded to update the property in question. A
 * instance of container class PropertyUpdater is created to hold
 * pointers to a set of updaters (instances of subclasses of
 * Updater). The conditions under which the properties must be updated
 * must be the same for all updaters pointed to by the PropertyUpdater
 * instance. For example, one instance of PropertyUpdater might handle
 * all properties that depent only on temperature, and another
 * properties that depend on composition.
 *
 * When updaters are added to a PropertyUpdater instance, an integer
 * is returned that can be used to refer to that updater.  The first
 * updater (number 0) is always a null updater that does nothing.
 *
 * PropertyUpdater acts as a 'switch'. When its method 'update(n)' is
 * called, the updater pointed to by the nth pointer in its internal
 * pointer array is invoked. The method 'need_update()' sets all
 * pointers in this array to point to the appropriate installed
 * updaters. When 'update(n)' is called, it invokes whatever updater
 * is pointed to, then sets the pointer to the null
 * updater. Subsequent calls to 'update(n)' will do nothing, until
 * 'need_update()' is called again. Note that 'update(n)' only sets
 * the nth pointer to the null updater; the rest continue to function
 * until 'update' is called with their index number. 'need_update',
 * however, resets all pointers to the non-null updaters.
 */

namespace Cantera {

    /**
     * Base class for updaters. This also serves as the null
     * updater. Specific updaters should be derived from this
     * class. Here's an example of a template for an updater that
     * calls method 'recompute()' of whatever object it is
     * initialized with:
     * @code
     * template<class S>
     * class UpdateMyProperty : public Updater {
     * public:
     *    UpdateMyProperty(S& s) :  m_s(s) {}
     *    virtual void update() { m_s.recompute(); }
     * private:
     *    S& m_s;
     * };
     * @code
     *
     * @ingroup updategroup
     */
    struct Updater {
        Updater() : count (0) {}
        virtual void update() {}
        long count;
    };


    /**
     * Property updater.
     * @ingroup updategroup
     */
    class PropertyUpdater {

    public:

        /**
         * Construct a new instance and install a null updater.
         */        
        PropertyUpdater() { 
            m_updaters.push_back( new Updater() );
            m_switch.push_back(0);
            m_number = 1;
        }

        
        /// Destructor. Does nothing.
        virtual ~PropertyUpdater() {}

        
        /**
         * Install an updater.
         * @param u pointer to an instance of a class derived from Updater.
         */
        int install(Updater* u) {
            m_updaters.push_back(u);
            m_switch.push_back(m_number);
            m_number++;
            need_update();
            return m_number - 1;
        }


        /**
         * Signal that an update is needed.
         */
        void need_update() {
            // reset all switches to the installed updaters
            for(int i=1; i < m_number; i++) m_switch[i] = i;
        }


        /**
         * Invoke the updater at position n. Depending on whether
         * need_update() or update(n) was called last, this will be
         * either the nth installed updater, or the null updater.
         */
        void update(int n) {
            m_updaters[m_switch[n]]->update();
            m_switch[n] = 0;                    // switch to the null updater
        }


        /**
         * Force an update of all properties, whether needed or not.
         */
        void force_update() {
                for (int n=1; n<m_number; n++) {
                        m_updaters[n]->update();
                }
        }
        
    private:

        vector<Updater*> m_updaters;
        vector_int   m_switch;
        int m_number;
    };

}

#endif






