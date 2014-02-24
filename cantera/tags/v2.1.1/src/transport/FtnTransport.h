/**
 *  @file FtnTransport.h
 *
 *  Customizable Fortran transport manager. This manager calls
 *  external Fortran functions to evaluate the transport
 *  properties. This is designed to be used to build custom transport
 *  managers in Fortran.
 */

// Copyright 2003  California Institute of Technology


#ifndef CT_FTNTRANSPORT_H
#define CT_FTNTRANSPORT_H

#include "cantera/transport/TransportBase.h"


/**
 * Change these definitions to change the names of the Fortran
 * procedures. If you want to define more than one custom Fortran
 * transport manager, copy this file, rename the class, and change
 * these definitions.
 *
 * The Fortran procedure names must follow the conventions of the Fortran
 * compiler. On most unix systems, this means the names must be lowercase,
 * and must include a trailing underscore.
 *
 */
#define __VISC__           visc_
#define __BULKVISC__       bvisc_
#define __TCON__           tcon_
#define __TDIFF__          tdiff_
#define __MULTIDIFF__      multidiff_
#define __MIXDIFF__        mixdiff_
#define __SIGMA__          sigma_
#define __GETMOBILITIES__  getmobilities_


extern "C" {

    doublereal __VISC__(doublereal* t, doublereal* p, doublereal* x);
    doublereal __BULKVISC__(doublereal* t, doublereal* p, doublereal* x);
    doublereal __TCON__(doublereal* t, doublereal* p, doublereal* x);

    void __TDIFF__(doublereal* t, doublereal* p, doublereal* x, doublereal* dt);
    void __MULTIDIFF__(doublereal* t, doublereal* p, doublereal* x,
                       integer* ld, doublereal* d);
    void __MIXDIFF__(doublereal* t, doublereal* p, doublereal* x, doublereal* d);
    void __BINDIFF__(doublereal* t, doublereal* p, doublereal* x, integer* ld, doublereal* d);
    doublereal __SIGMA__(doublereal* t, doublereal* p, doublereal* x);

    doublereal __GETMOBILITIES__(doublereal* t, doublereal* p,
                                 doublereal* x, doublereal* mobil);

}

namespace Cantera
{

/**
 * A class that calls external Fortran functions to evaluate
 * transport properties.
 * @deprecated Broken and unused
 */
class FtnTransport : public Transport
{

public:

    FtnTransport(int model, thermo_t* thermo) : Transport(thermo) {
        warn_deprecated("FtnTransport", "This class will be removed in Cantera 2.2.");
        m_model = model;
        m_x.resize(m_thermo->nSpecies(), 0.0);
        updateTPX();
    }

    virtual int model() {
        return cFtnTransport + m_model;
    }

    virtual doublereal viscosity() {
        updateTPX();
        return __VISC__(&m_temp, &m_pres, m_x.begin());
    }

    virtual doublereal bulkViscosity() {
        updateTPX();
        return __BULKVISC__(&m_temp, &m_pres, m_x.begin());
    }

    virtual doublereal thermalConductivity() {
        updateTPX();
        return __TCON__(&m_temp, &m_pres, m_x.begin());
    }

    virtual doublereal electricalConductivity() {
        updateTPX();
        return __SIGMA__(&m_temp, &m_pres, m_x.begin());
    }

    virtual void getMobilities(doublereal* mobil) {
        updateTPX();
        __GETMOBILITIES__(&m_temp, &m_pres, m_x.begin(), mobil);
    }


    virtual void getThermalDiffCoeffs(doublereal* dt) {
        updateTPX();
        __TDIFF__(&m_temp, &m_pres, m_x.begin(), dt);
    }

    virtual void getBinaryDiffCoeffs(int ld, doublereal* d) {
        updateTPX();
        integer ldd = ld;
        __BINDIFF__(&m_temp, &m_pres, m_x.begin(), &ldd, d);
    }

    virtual void getMultiDiffCoeffs(int ld, doublereal* d) {
        updateTPX();
        integer ldd = ld;
        __MULTIDIFF__(&m_temp, &m_pres, m_x.begin(), &ldd, d);
    }

    virtual void getMixDiffCoeffs(doublereal* d) {
        updateTPX();
        __MIXDIFF__(&m_temp, &m_pres, m_x.begin(), d);
    }


private:

    void updateTPX() {
        m_temp = m_thermo->temperature();
        m_pres = m_thermo->pressure();
        m_thermo->getMoleFractions(m_x.begin());
    }
    doublereal m_temp;
    doublereal m_pres;
    vector_fp  m_x;
    int m_model;

};

}
#endif






