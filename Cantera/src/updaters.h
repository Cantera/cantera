#ifndef CT_UPDATERS_H
#define CT_UPDATERS_H

#include "PropertyUpdater.h"

namespace Cantera {

    //--------------------------------------------------------
    //         Property Updaters
    //--------------------------------------------------------

    /**
     * Invokes method 'update_T' of the object it is initialized with.
     */
    template<class S>
    struct T_Updater : public Updater {
        T_Updater(S& s) : m_s(s) {}
        void update() { m_s.update_T(); }
        S& m_s;
    };


    /**
     * Invokes method 'updateMoleFractions' of the object it is
     * initialized with.
     */
    template<class S>
    struct UpdateMoleFractions : public Updater {
        UpdateMoleFractions(S& s) : Updater(), m_s(s) {}
        void update() { m_s.updateMoleFractions(); }
        S& m_s;
    };

    /**
     * Invokes method 'updateMW' of the object it is
     * initialized with.
     */
    template<class S>
    struct UpdateMolWt : public Updater {
        UpdateMolWt(S& s) : Updater(), m_s(s) {}
        void update() { m_s.updateMW(); }
        S& m_s;
    };

    /**
     * Updater responsible for updating the species standard-state
     * thermodynamic properties.
     * @ingroup updategroup
     */
    template<class S>
    struct UpdateThermo : public Updater {
        UpdateThermo(S& s) : Updater(), m_s(s) {}
        void update() { m_s._updateThermo(); }

        S& m_s;
    };

    /**
     * Updater responsible for updating the temperature-dependent
     * parts of the transport properties.
     * @ingroup updategroup
     */
    template<class S>
    struct UpdateTransport_T : public Updater {
        UpdateTransport_T(S& s) : Updater(), m_s(s) {}
        void update() { m_s._update_transport_T(); }
        S& m_s;
    };


    /**
     * Updater responsible for updating the concentration-dependent
     * parts of the transport properties.
     * @ingroup updategroup
     */
    template<class S>
    struct UpdateTransport_C : public Updater {
        UpdateTransport_C(S& s) : Updater(), m_s(s) {}
        void update() { m_s._update_transport_C(); }
        S& m_s;
    };    


    /**
     * Template for an updater subclass that calls method _updateRates_T()
     * of the object it is initialized with.
     * @ingroup updategroup   
     */
    template<class S>
    struct UpdateRates_T : public Updater {
        UpdateRates_T(S& s) : Updater(), m_s(s) {}
        void update() { m_s._update_rates_T(); }
        S& m_s;
    };

    /**
     * Template for an updater subclass that calls method _updateRates_C()
     * of the object it is initialized with.
     * @ingroup updategroup  
     */
    template<class S>
    struct UpdateRates_C : public Updater {
        UpdateRates_C(S& s) : Updater(), m_s(s) {}
        void update() { m_s._update_rates_C(); }
        S& m_s;
    };


}

#endif


