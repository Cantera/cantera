#ifndef CT_MIX_UTILS_H
#define CT_MIX_UTILS_H

#include "ctexceptions.h"

namespace Cantera {

    doublereal quadInterp(doublereal x0, doublereal* x, doublereal* y);

//     /** Set the temperature (K), pressure (Pa), and mole fractions.  */
//     template<class S>
//     void set_TPX(S& s, doublereal t, doublereal p, doublereal* x) {
//         s.setMoleFractions(x); s.setTemperature(t); s.setPressure(p);
//     }

//     template<class S>
//     void set_HP(S& s, doublereal h, doublereal p, doublereal tol) {
//         doublereal dt;
//         s.setPressure(p);
//         for (int n = 0; n < 20; n++) {
//             dt = (h - s.enthalpy_mass())/s.cp_mass();
//             if (dt > 100.0) dt = 100.0;
//             else if (dt < -100.0) dt = -100.0; 
//             s.setState_TP(s._temp() + dt, p);
//             if (fabs(dt) < tol) {
//                 return;
//             }
//         }
//         throw CanteraError("set_HP","no convergence. dt = " + fp2str(dt));
//     }

//     template<class S>
//     void set_UV(S& s, doublereal u, doublereal v, 
// 				    doublereal tol) {
//         doublereal dt;
//         s.setDensity(1.0/v);
//         for (int n = 0; n < 20; n++) {
//             dt = (u - s.intEnergy_mass())/s.cv_mass();
//             if (dt > 100.0) dt = 100.0;
//             else if (dt < -100.0) dt = -100.0; 
//             s.setTemperature(s._temp() + dt);
//             if (fabs(dt) < tol) return;
//         }
//         throw CanteraError("set_UV","no convergence. dt = " + fp2str(dt));
//     }

//     template<class S>
//     void set_SP(S& s, doublereal entropy, doublereal p, 
//         doublereal tol) {
//         doublereal dt;
//         s.setPressure(p);
//         for (int n = 0; n < 20; n++) {
//             dt = (entropy - s.entropy_mass())*s.temperature()/s.cp_mass();
//             if (dt > 100.0) dt = 100.0;
//             else if (dt < -100.0) dt = -100.0; 
//             s.setState_TP(s._temp() + dt, p);
//             if (fabs(dt) < tol) return;
//         }
//         throw CanteraError("set_SP","no convergence. dt = " + fp2str(dt));
//     }

//     template<class S>
//     void set_SV(S& s, doublereal entropy, doublereal v, 
//         doublereal tol) {
//         doublereal dt;
//         s.setDensity(1.0/v);
//         for (int n = 0; n < 20; n++) {
//             dt = (entropy - s.entropy_mass())*s.temperature()/s.cv_mass();
//             if (dt > 100.0) dt = 100.0;
//             else if (dt < -100.0) dt = -100.0; 
//             s.setTemperature(s._temp() + dt);
//             if (fabs(dt) < tol) return;
//         }
//         throw CanteraError("set_SV","no convergence. dt = " + fp2str(dt));
//     }

    template<class S1, class S2>
    void mapSpeciesData(const S1& s1, const S2& s2, const doublereal* data1,
        doublereal* data2) {
        int n1 = s1.nSpecies();
        int n2 = s2.nSpecies();
        int n, m2;

        // zero out the destination array
        for (n = 0; n < n2; n++) data2[n] = 0.0;

        // copy 
        for (n = 0; n < n1; n++) {
            m2 = s2.speciesIndex(s1.speciesName(n));
            if (m2 >= 0) data2[m2] = data1[n];
        }
    }
}

#endif
