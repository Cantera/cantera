/**
 * @file fitPoly.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_FITPOLY_H
#define CT_FITPOLY_H

#include "../Cantera/src/polyfit.h"

namespace Cantera {

    inline doublereal poly_enthalpy_RT(doublereal temp, 
        int n, doublereal* c, doublereal h0) {
        doublereal tpwr = 1.0;
        int j;
        doublereal h = 0.0;
        for (j = 0; j <= n; j++) {
            h += c[j]*tpwr/(j+1);
            tpwr *= temp;
        }
        return h + h0/temp;
    }

    inline doublereal poly_entropy_R(doublereal temp, 
        int n, doublereal* c, doublereal s0) {
        doublereal tpwr = temp;
        int j;
        doublereal s = c[0]*log(temp);
        for (j = 1; j <= n; j++) {
            s += c[j]*tpwr/(j);
            tpwr *= temp;
        }
        return s + s0;
    }
    
    template<class M>
    void fitPolynomial(int n, M& mix, int fittype, vector<vector_fp>& coeffs, 
        doublereal tmin = -1.0, doublereal tmax = -1.0) {
        int i, k;
        doublereal temp;
        if (tmin < 0.0) tmin = mix.minTemp();
        if (tmax < 0.0) tmax = mix.maxTemp();

        // generate data
        int np = int((tmax - tmin)/20.0 + 1);
        if (np < 2*n+2) np = 2*n+2;
        doublereal dt = (tmax - tmin)/(np - 1);
        
        int nsp = mix.nSpecies();

        coeffs.resize(nsp);

        vector<vector_fp> cp(nsp);
        vector_fp x(np), w(np);
        for (k = 0; k < nsp; k++) {
            cp[k].resize(np);
        }

        mix.setTemperature(298.15);
        vector_fp h0 = mix.enthalpy_RT();
        vector_fp s0 = mix.entropy_R();
        
        for (i = 0; i < np; i++) {
            temp = tmin + i*dt;
            mix.setTemperature(temp);
            const vector_fp& cpr = mix.cp_R();
            switch (fittype) {
            case 0: x[i] = temp; break;
            case 1: x[i] = 1.0/temp; break;
            case 2: x[i] = log(temp); break;
            default: 
                cout << " unknown fit type (" << fittype << ")" << endl;
                return;
            }
            w[i] = -1.0;
            for (k = 0; k < nsp; k++) {
                cp[k][i] = cpr[k];
            }
        }
        doublereal err;
        for (k = 0; k < nsp; k++) {
            coeffs[k].resize(n + 3);
            err = polyfit(np, x.begin(), cp[k].begin(), w.begin(),
                n, n, 0.0, coeffs[k].begin()+2);
            coeffs[k][0] = 298.15 * (h0[k] - 
                poly_enthalpy_RT(298.15, n, coeffs[k].begin()+2, 0.0));
            coeffs[k][1] = (s0[k] 
                - poly_entropy_R(298.15, n, coeffs[k].begin()+2, 0.0));
        }
    }
}

#undef TESTIT
#ifdef TESTIT

#include "Cantera.h"
#include "../Cantera/src/IdealGasMix.h"

main() {
    IdealGasMix mix("gri30.inp");
    int fittype = 0;
    int n;
    cout << " enter n: ";
    cin >> n;
    fitPolynomial(n, mix, fittype); 
}

#endif

#endif
