#ifndef CT_UNITS_H
#define CT_UNITS_H

#include "ct_defs.h"

namespace Cantera {

    class Unit {
    public:

        static Unit* units() {
            if (!__u) __u = new Unit;
            return __u;
        }

        virtual ~Unit() {
            delete __u;
            __u = 0;
        }

        doublereal toSI(string units) {
            if (units == "") return 1.0;
            doublereal f = 1.0, fctr;
            int tsize;
            string u = units, tok, tsub;
            int k;
            char action = '-';
            //if (units[0] == '/') action = '/';
            while (1 > 0) {
                k = u.find_first_of("/-");
                if (k >= 0)
                    tok = u.substr(0,k);
                else
                    tok = u;
                tsize = tok.size();
                if (tsize == 0) 
                    fctr = 1.0;
                else if (tok[tsize - 1] == '2') {
                    tsub = tok.substr(0,tsize-1);
                    fctr = m_u[tsub];
                    fctr *= fctr;
                }
                else if (tok[tsize - 1] == '3') {
                    tsub = tok.substr(0,tsize-1);
                    fctr = m_u[tsub];
                    fctr *= fctr*fctr;
                }
                else if (tok[tsize - 1] == '4') {
                    tsub = tok.substr(0,tsize-1);
                    fctr = m_u[tsub];
                    fctr *= fctr*fctr*fctr;
                }
                else if (tok[tsize - 1] == '5') {
                    tsub = tok.substr(0,tsize-1);
                    fctr = m_u[tsub];
                    fctr *= fctr*fctr*fctr*fctr;
                }
                else if (tok[tsize - 1] == '6') {
                    tsub = tok.substr(0,tsize-1);
                    fctr = m_u[tsub];
                    fctr *= fctr*fctr*fctr*fctr*fctr;
                }
                else {
                    tsub = tok;
                    fctr = m_u[tok];
                }

                if (fctr == 0) 
                    throw CanteraError("toSI","unknown unit: "+tsub);
                if (action == '-') f *= fctr;
                else if (action == '/') f /= fctr;
                if (k < 0) break;
                action = u[k];
                u = u.substr(k+1,u.size());
            }
            return f;
        }

    private:

        static Unit* __u;
        map<string, doublereal> m_u;
        Unit(){

            // length
            m_u["m"]    = 1.0;
            m_u["cm"]   = 0.01;
            m_u["km"]   = 1.0e3;
            m_u["mm"]   = 1.0e-3;
            m_u["micron"]   = 1.0e-6;
            m_u["A"]   = 1.0e-10;

            // energy
            m_u["J"]        = 1.0;
            m_u["kJ"]       = 1.0e3;
            m_u["cal"]      = 4.184;
            m_u["kcal"]     = 4184.0;

            // quantity
            m_u["mol"]      = 1.0e-3;
            m_u["mole"]      = 1.0e-3;
            m_u["kmol"]     = 1.0;
            m_u["molec"]    = 1.0/Avogadro;

            // temperature
            m_u["K"]        = 1.0;
            m_u["C"]        = 1.0;

            // mass
            m_u["g"]        = 1.0e-3;
            m_u["kg"]       = 1.0;

            // pressure
            m_u["atm"]      = 1.01325e5;
            m_u["bar"]      = 1.0e5;
            m_u["Pa"]       = 1.0;

            // time
            m_u["s"]        = 1.0;
            m_u["min"]      = 60.0;
            m_u["hr"]       = 3600.0;
            m_u["ms"]       = 0.001;
        }
    };
}

#endif
