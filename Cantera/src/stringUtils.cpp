/**
 *  @file stringUtils.cpp
 *
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
//#include "Phase.h"
#include "ThermoPhase.h"

//#include "IdealGasMix.h"

#include <stdio.h>

namespace Cantera {


    /**
     * Convert a floating point number to a string using sprintf.
     */
    string fp2str(double x, string fmt) {
        char buf[30];
        sprintf(buf, fmt.c_str(), x);
        return string(buf);
    }

    /**
     * Convert an integer number to a string using sprintf.
     */
    string int2str(int n, string fmt) {
        char buf[30];
        sprintf(buf, fmt.c_str(), n);
        return string(buf);
    }

    string lowercase(string s) {
        int n = s.size();
        string lc(s);
        for (int i = 0; i < n; i++) lc[i] = tolower(s[i]);
        return lc;
    }

    /** 
     * Strip white space.
     */
    string stripws(string s) {
        int i;
        bool sempty = true;
        int n = s.size();
        string ss = "";
        for (i = 0; i < n; i++) {
            if (s[i] != ' ' && isprint(s[i])) {
                ss += s[i];
                sempty = false;
            }
            else if (!sempty) break;
        }
        return ss;
    }

    /** 
     * Strip non-printing characters.
     */
    string stripnonprint(string s) {
        int i;
        int n = s.size();
        string ss = "";
        for (i = 0; i < n; i++) {
            if (isprint(s[i])) {
                ss += s[i];
            }
        }
        return ss;
    }

            
    /**
     * Parse a composition string.
     */
    void parseCompString(const string ss, compositionMap& x) {
        string s = ss;
        int icolon, ibegin, iend;
        string name, num, nm;
        do {
            ibegin = s.find_first_not_of(", ;\n\t");
            if (ibegin >= 0) {
                s = s.substr(ibegin,s.size());
                icolon = s.find(':');
                iend = s.find_first_of(", ;\n\t");
                //icomma = s.find(',');
                if (icolon > 0) {
                    name = s.substr(0, icolon);
                    if (iend > 0) {
                        num = s.substr(icolon+1, iend-icolon);
                        s = s.substr(iend+1, s.size());
                    }
                    else {
                        num = s.substr(icolon+1, s.size());
                        s = "";
                    }
                    nm = stripws(name);
                    if (x[nm] == 0.0) throw CanteraError("parseCompString",
                        "unknown species "+nm);
                    x[nm] = atof(num.c_str());
                }
                else s = "";
            }
        }
        while (s != "");
    }

    int fillArrayFromString(const string& str, doublereal* a, char delim) {
        int iloc;
        int count = 0;
        string num, s;
        s = str;
        while (s.size() > 0) {
            iloc = s.find(delim);
            if (iloc > 0) {
                num = s.substr(0, iloc);
                s = s.substr(iloc+1,s.size());
            }
            else {
                num = s;
                s = "";
            }
            a[count] = atof(num.c_str());
            count++;
        }
        return count;
    }

    /**
     * Format a summary of the mixture state for output.
     */           
    string report(const ThermoPhase& th, bool show_thermo) {

        char p[200];
        string s = "";

        sprintf(p, "\n       temperature    %12.6g  K\n", th.temperature());
        s += p;
        sprintf(p, "          pressure    %12.6g  Pa\n", th.pressure());
        s += p;
        sprintf(p, "           density    %12.6g  kg/m^3\n", th.density());
        s += p;
        sprintf(p, "  mean mol. weight    %12.6g  amu\n", th.meanMolecularWeight());
        s += p;

        if (show_thermo) {
        sprintf(p, "\n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n", 
            th.enthalpy_mass(), th.enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n", 
            th.intEnergy_mass(), th.intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n", 
            th.entropy_mass(), th.entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n", 
            th.gibbs_mass(), th.gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n", 
            th.cp_mass(), th.cp_mole());
        s += p;
        sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
            th.cv_mass(), th.cv_mole());
        s += p;
        }

        int kk = th.nSpecies();
        array_fp x(kk);
        array_fp y(kk);
        th.getMoleFractions(x.begin());
        th.getMassFractions(y.begin());

        int k;

        sprintf(p, "\n                           X                 Y   \n");
        s += p;
        sprintf(p, "                     -------------     ------------\n");
        s += p;
        for (k = 0; k < kk; k++) {
                sprintf(p, "%18s   %12.6e     %12.6e\n", 
                    th.speciesName(k).c_str(), x[k], y[k]);
                s += p;
        }
        return s;
    }

    /**
     * Format a composition list for output.
     */           
    string formatCompList(const Phase& mix, int xyc) {

        const doublereal Threshold = 1.e-20;

        char p[200];
        string s = "";
        int kk = mix.nSpecies();
        array_fp zz(kk);
        switch (xyc) {
        case 0: mix.getMoleFractions(zz.begin()); break;
        case 1: mix.getMassFractions(zz.begin()); break;
        case 2: mix.getConcentrations(zz.begin()); break;
        default: return "error: xyc must be 0, 1, or 2";
        }

        doublereal z;
        int k;
        for (k = 0; k < kk; k++) {
            z = fabs(zz[k]);
            if (z < Threshold) zz[k] = 0.0;
        }

        for (k = 0; k < kk; k++) {
            sprintf(p, "%18s\t %12.6e\n", mix.speciesName(k).c_str(), 
                zz[k]);
            s += p;
        }
        return s;
    }


    /**
     *  Get the file name without the path or extension
     */
    string getFileName(const string& path) {
        string file;
        size_t idot = path.find_last_of('.');
        size_t islash = path.find_last_of('/');
        if (idot > 0 && idot < path.size()) {
            if (islash > 0 && islash < idot) {
                file = path.substr(islash+1, idot-islash-1);
            }
            else {
                file = path.substr(0,idot);
            }
        }
        else {
            file = path;
        }       
        return file;
    }


    /**
     *  Generate a logfile name based on an input file name
     */
    string logfileName(const string& infile) {
        string logfile;
        size_t idot = infile.find_last_of('.');
        size_t islash = infile.find_last_of('/');
        if (idot > 0 && idot < infile.size()) {
            if (islash > 0 && islash < idot) {
                logfile = infile.substr(islash+1, idot-islash-1) + ".log";
            }
            else {
                logfile = infile.substr(0,idot) + ".log";
            }
        }
        else {
            logfile = infile + ".log";
        }       
        return logfile;
    }
}
