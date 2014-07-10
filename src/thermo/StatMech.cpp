/**
 *  @file StatMech.cpp
 *  \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink
 */

// Copyright 2007  Sandia National Laboratories

#include "cantera/thermo/StatMech.h"
#include <iostream>

namespace Cantera
{
StatMech::StatMech() {
    warn_deprecated("class StatMech", "To be removed after Cantera 2.2");
}

StatMech::StatMech(int n, doublereal tlow, doublereal thigh,
                   doublereal pref,
                   const doublereal* coeffs,
                   const std::string& my_name) :
    SpeciesThermoInterpType(n, tlow, thigh, pref),
    sp_name(my_name)
{
    // should error on zero -- cannot take ln(0)
    if (m_lowT <= 0.0) {
        throw CanteraError("Error in StatMech.cpp",
                           " Cannot take 0 tmin as input. \n\n");
    }
    buildmap();
}

StatMech::StatMech(const StatMech& b) :
    SpeciesThermoInterpType(b)
{
}

StatMech& StatMech::operator=(const StatMech& b)
{
    if (&b != this) {
        SpeciesThermoInterpType::operator=(b);
    }
    return *this;
}

SpeciesThermoInterpType*
StatMech::duplMyselfAsSpeciesThermoInterpType() const
{
    return new StatMech(*this);
}

int StatMech::reportType() const
{
    return STAT;
}

int StatMech::buildmap()
{

    // build vector of strings
    std::vector<std::string> SS;

    // now just iterate over name map to place each
    // string in a key

    SS.push_back("Air");
    SS.push_back("CPAir");
    SS.push_back("Ar");
    SS.push_back("Ar+");
    SS.push_back("C");
    SS.push_back("C+");
    SS.push_back("C2");
    SS.push_back("C2H");
    SS.push_back("C2H2");
    SS.push_back("C3");
    SS.push_back("CF");
    SS.push_back("CF2");
    SS.push_back("CF3");
    SS.push_back("CF4");
    SS.push_back("CH");
    SS.push_back("CH2");
    SS.push_back("CH3");
    SS.push_back("CH4");
    SS.push_back("Cl");
    SS.push_back("Cl2");
    SS.push_back("CN");
    SS.push_back("CN+");
    SS.push_back("CO");
    SS.push_back("CO+");
    SS.push_back("CO2");
    SS.push_back("F");
    SS.push_back("F2");
    SS.push_back("H");
    SS.push_back("H+");
    SS.push_back("H2");
    SS.push_back("H2+");
    SS.push_back("H2O");
    SS.push_back("HCl");
    SS.push_back("HCN");
    SS.push_back("He");
    SS.push_back("He+");
    SS.push_back("N");
    SS.push_back("N+");
    SS.push_back("N2");
    SS.push_back("CPN2");
    SS.push_back("N2+");
    SS.push_back("Ne");
    SS.push_back("NCO");
    SS.push_back("NH");
    SS.push_back("NH+");
    SS.push_back("NH2");
    SS.push_back("NH3");
    SS.push_back("NO");
    SS.push_back("NO+");
    SS.push_back("NO2");
    SS.push_back("O");
    SS.push_back("O+");
    SS.push_back("O2");
    SS.push_back("O2+");
    SS.push_back("OH");
    SS.push_back("Si");
    SS.push_back("SiO");
    SS.push_back("e");

    // now place each species in a map
    size_t ii;
    for (ii=0; ii < SS.size(); ii++) {
        name_map[SS[ii]]=(new species);

        // init to crazy defaults
        name_map[SS[ii]]->nvib       = -1;
        name_map[SS[ii]]->cfs        = -1;
        name_map[SS[ii]]->mol_weight = -1;

        name_map[SS[ii]]->theta[0]   =0.0;
        name_map[SS[ii]]->theta[1]   =0.0;
        name_map[SS[ii]]->theta[2]   =0.0;
        name_map[SS[ii]]->theta[3]   =0.0;
        name_map[SS[ii]]->theta[4]   =0.0;
    }

    // now set all species information

    // build Air
    name_map["Air"]->cfs = 2.5;
    name_map["Air"]->mol_weight=28.96;
    name_map["Air"]->nvib=0;

    // build CPAir
    name_map["CPAir"]->cfs = 2.5;
    name_map["CPAir"]->mol_weight=28.96;
    name_map["CPAir"]->nvib=0;

    // build Ar
    name_map["Ar"]->cfs = 1.5;
    name_map["Ar"]->mol_weight=39.944;
    name_map["Ar"]->nvib=0;

    // build Ar+
    name_map["Ar+"]->cfs = 1.5;
    name_map["Ar+"]->mol_weight=39.94345;
    name_map["Ar+"]->nvib=0;

    // build C
    name_map["C"]->cfs = 1.5;
    name_map["C"]->mol_weight=12.011;
    name_map["C"]->nvib=0;

    // build C+
    name_map["C+"]->cfs = 1.5;
    name_map["C+"]->mol_weight=12.01045;
    name_map["C+"]->nvib=0;

    // C2
    name_map["C2"]->cfs=2.5;
    name_map["C2"]->mol_weight=24.022;
    name_map["C2"]->nvib=1;
    name_map["C2"]->theta[0]=2.6687e3;

    // C2H
    name_map["C2H"]->cfs=2.5;
    name_map["C2H"]->mol_weight=25.03;
    name_map["C2H"]->nvib=3;
    name_map["C2H"]->theta[0]=5.20100e+03;
    name_map["C2H"]->theta[1]=7.20000e+03;
    name_map["C2H"]->theta[2]=2.66100e+03;

    // C2H2
    name_map["C2H2"]->cfs=2.5;
    name_map["C2H2"]->mol_weight=26.038;
    name_map["C2H2"]->nvib=5;
    name_map["C2H2"]->theta[0]=4.85290e+03;
    name_map["C2H2"]->theta[1]=2.84000e+03;
    name_map["C2H2"]->theta[2]=4.72490e+03;
    name_map["C2H2"]->theta[3]=8.81830e+02;
    name_map["C2H2"]->theta[4]=1.05080e+03;

    // C3
    name_map["C3"]->cfs=2.5;
    name_map["C3"]->mol_weight=36.033;
    name_map["C3"]->nvib=3;
    name_map["C3"]->theta[0]=1.84500e+03;
    name_map["C3"]->theta[1]=7.78700e+02;
    name_map["C3"]->theta[2]=3.11760e+03;

    // CF
    name_map["CF"]->cfs=2.5;
    name_map["CF"]->mol_weight=31.00940;
    name_map["CF"]->nvib=1;
    name_map["CF"]->theta[0]=1.88214e+03;

    // CF2
    name_map["CF2"]->cfs=3;
    name_map["CF2"]->mol_weight=50.00780;
    name_map["CF2"]->nvib=3;
    name_map["CF2"]->theta[0]=1.76120e+03;
    name_map["CF2"]->theta[1]=9.56820e+02;
    name_map["CF2"]->theta[2]=1.60000e+03;

    // CF3
    name_map["CF3"]->cfs=3;
    name_map["CF3"]->mol_weight=69.00620;
    name_map["CF3"]->nvib=4;
    name_map["CF3"]->theta[0]=1.56800e+03;
    name_map["CF3"]->theta[1]=1.00900e+03;
    name_map["CF3"]->theta[2]=1.81150e+03;
    name_map["CF3"]->theta[3]=7.36680e+02;

    // CF4
    name_map["CF4"]->cfs=3;
    name_map["CF4"]->mol_weight=88.00460;
    name_map["CF4"]->nvib=4;
    name_map["CF4"]->theta[0]=1.30720e+03;
    name_map["CF4"]->theta[1]=6.25892e+02;
    name_map["CF4"]->theta[2]=1.84540e+03;
    name_map["CF4"]->theta[3]=9.08950e+02;

    // CH
    name_map["CH"]->cfs=2.5;
    name_map["CH"]->mol_weight=13.01900;
    name_map["CH"]->nvib=1;
    name_map["CH"]->theta[0]=4.11290e+03;

    // CH2
    name_map["CH2"]->cfs=3;
    name_map["CH2"]->mol_weight=14.02700;
    name_map["CH2"]->nvib=3;
    name_map["CH2"]->theta[0]=4.31650e+03;
    name_map["CH2"]->theta[1]=1.95972e+03;
    name_map["CH2"]->theta[2]=4.60432e+03;

    // CH3
    name_map["CH3"]->cfs=3;
    name_map["CH3"]->mol_weight=15.03500;
    name_map["CH3"]->nvib=4;
    name_map["CH3"]->theta[0]=4.31650e+03;
    name_map["CH3"]->theta[1]=8.73370e+02;
    name_map["CH3"]->theta[2]=4.54960e+03;
    name_map["CH3"]->theta[3]=2.01150e+03;

    // CH4
    name_map["CH4"]->cfs=3;
    name_map["CH4"]->mol_weight=16.04300;
    name_map["CH4"]->nvib=4;
    name_map["CH4"]->theta[0]=4.19660e+03;
    name_map["CH4"]->theta[1]=2.20620e+03;
    name_map["CH4"]->theta[2]=4.34450e+03;
    name_map["CH4"]->theta[3]=1.88600e+03;

    // Cl
    name_map["Cl"]->cfs=1.5;
    name_map["Cl"]->mol_weight=35.45300;
    name_map["Cl"]->nvib=0;

    // Cl2
    name_map["Cl2"]->cfs=2.5;
    name_map["Cl2"]->mol_weight=70.96;
    name_map["Cl2"]->nvib=1;
    name_map["Cl2"]->theta[0]=8.05355e+02;

    // CN
    name_map["CN"]->cfs=2.5;
    name_map["CN"]->mol_weight=26.01900;
    name_map["CN"]->nvib=1;
    name_map["CN"]->theta[0]=2.97610e+03;

    // CN+
    name_map["CN+"]->cfs=2.5;
    name_map["CN+"]->mol_weight=26.01845;
    name_map["CN+"]->nvib=1;
    name_map["CN+"]->theta[0]=2.92520e+03;

    // CO
    name_map["CO"]->cfs=2.5;
    name_map["CO"]->mol_weight=28.01100;
    name_map["CO"]->nvib=1;
    name_map["CO"]->theta[0]=3.12200e+03;

    // CO+
    name_map["CO+"]->cfs=2.5;
    name_map["CO+"]->mol_weight=28.01045;
    name_map["CO+"]->nvib=1;
    name_map["CO+"]->theta[0]=3.18800e+03;

    // CO2
    name_map["CO2"]->cfs=2.5;
    name_map["CO2"]->mol_weight=44.01100;
    name_map["CO2"]->nvib=3;
    name_map["CO2"]->theta[0]=1.91870e+03;
    name_map["CO2"]->theta[1]=9.59660e+02;
    name_map["CO2"]->theta[2]=3.38210e+03;

    // F
    name_map["F"]->cfs=1.5;
    name_map["F"]->mol_weight=18.99840;
    name_map["F"]->nvib=0;

    // F2
    name_map["F2"]->cfs=2.5;
    name_map["F2"]->mol_weight=37.99680;
    name_map["F2"]->nvib=1;
    name_map["F2"]->theta[0]=1.32020e+03;

    // H
    name_map["H"]->cfs=1.5;
    name_map["H"]->mol_weight=1;
    name_map["H"]->nvib=0;

    // H+
    name_map["H+"]->cfs=1.5;
    name_map["H+"]->mol_weight=1.00745;
    name_map["H+"]->nvib=0;

    // H2
    name_map["H2"]->cfs=2.5;
    name_map["H2"]->mol_weight=2.01600;
    name_map["H2"]->nvib=1;
    name_map["H2"]->theta[0]=6.33140e+03;

    // H2+
    name_map["H2+"]->cfs=2.5;
    name_map["H2+"]->mol_weight=2.01545;
    name_map["H2+"]->nvib=1;
    name_map["H2+"]->theta[0]=3.34280e+03;

    // H2O
    name_map["H2O"]->cfs=3.0;
    name_map["H2O"]->mol_weight=18.01600;
    name_map["H2O"]->nvib=3;
    name_map["H2O"]->theta[0]=5.26130e+03;
    name_map["H2O"]->theta[1]=2.29460e+03;
    name_map["H2O"]->theta[2]=5.40395e+03;

    // HCl
    name_map["HCl"]->cfs=2.5;
    name_map["HCl"]->mol_weight=36.46100;
    name_map["HCl"]->nvib=1;
    name_map["HCl"]->theta[0]=4.30330e+03;

    // HCN
    name_map["HCN"]->cfs=2.5;
    name_map["HCN"]->mol_weight=27.02700;
    name_map["HCN"]->nvib=3;
    name_map["HCN"]->theta[0]=3.01620e+03;
    name_map["HCN"]->theta[1]=1.02660e+03;
    name_map["HCN"]->theta[2]=4.76450e+03;

    // He
    name_map["He"]->cfs=1.5;
    name_map["He"]->mol_weight=4.00300;
    name_map["He"]->nvib=0;

    // He+
    name_map["He+"]->cfs=1.5;
    name_map["He+"]->mol_weight=4.00245;
    name_map["He+"]->nvib=0;

    // N
    name_map["N"]->cfs=1.5;
    name_map["N"]->mol_weight=14.008;
    name_map["N"]->nvib=0;

    // Ne
    name_map["Ne"]->cfs=1.5;
    name_map["Ne"]->mol_weight=20.17900;
    name_map["Ne"]->nvib=0;

    // N+
    name_map["N+"]->cfs=1.5;
    name_map["N+"]->mol_weight=14.00745;
    name_map["N+"]->nvib=0;

    // N2
    name_map["N2"]->cfs=2.5;
    name_map["N2"]->mol_weight=28.01600;
    name_map["N2"]->nvib=1;
    name_map["N2"]->theta[0]=3.39500e+03;

    // N2+
    name_map["N2+"]->cfs=2.5;
    name_map["N2+"]->mol_weight=28.01545;
    name_map["N2+"]->nvib=1;
    name_map["N2+"]->theta[0]=3.17580e+03;

    // CPN2
    name_map["CPN2"]->cfs=2.5;
    name_map["CPN2"]->mol_weight=28.01600;
    name_map["CPN2"]->nvib=0;

    // NCO
    name_map["NCO"]->cfs=2.5;
    name_map["NCO"]->mol_weight=42.01900;
    name_map["NCO"]->nvib=3;
    name_map["NCO"]->theta[0]=1.83600e+03;
    name_map["NCO"]->theta[1]=7.67100e+02;
    name_map["NCO"]->theta[2]=2.76800e+03;

    // NH
    name_map["NH"]->cfs=2.5;
    name_map["NH"]->mol_weight=15.01600;
    name_map["NH"]->nvib=1;
    name_map["NH"]->theta[0]=4.72240e+03;

    // NH+
    name_map["NH+"]->cfs=2.5;
    name_map["NH+"]->mol_weight=15.01545;
    name_map["NH+"]->nvib=0;

    // NH2
    name_map["NH2"]->cfs=2.5;
    name_map["NH2"]->mol_weight=16.02400;
    name_map["NH2"]->nvib=0;

    // NH3
    name_map["NH3"]->cfs=2.5;
    name_map["NH3"]->mol_weight=17.03200;
    name_map["NH3"]->nvib=4;
    name_map["NH3"]->theta[0]=4.78100e+03;
    name_map["NH3"]->theta[1]=1.47040e+03;
    name_map["NH3"]->theta[2]=4.95440e+03;
    name_map["NH3"]->theta[3]=2.34070e+03;

    // NO
    name_map["NO"]->cfs=2.5;
    name_map["NO"]->mol_weight=30.00800;
    name_map["NO"]->nvib=1;
    name_map["NO"]->theta[0]=2.81700e+03;

    // NO+
    name_map["NO+"]->cfs=2.5;
    name_map["NO+"]->mol_weight=30.00745;
    name_map["NO+"]->nvib=1;
    name_map["NO+"]->theta[0]=3.42100e+03;

    // NO2
    name_map["NO2"]->cfs=3;
    name_map["NO2"]->mol_weight=46.00800;
    name_map["NO2"]->nvib=3;
    name_map["NO2"]->theta[0]=1.07900e+03;
    name_map["NO2"]->theta[1]=1.90000e+03;
    name_map["NO2"]->theta[2]=2.32700e+03;

    // O
    name_map["O"]->cfs=1.5;
    name_map["O"]->mol_weight=16.000;
    name_map["O"]->nvib=0;

    // O+
    name_map["O+"]->cfs=1.5;
    name_map["O+"]->mol_weight=15.99945;
    name_map["O+"]->nvib=0;

    // O2
    name_map["O2"]->cfs=2.5;
    name_map["O2"]->mol_weight=32.00000;
    name_map["O2"]->nvib=1;
    name_map["O2"]->theta[0]=2.23900e+03;

    // O2
    name_map["O2+"]->cfs=2.5;
    name_map["O2+"]->mol_weight=31.99945;
    name_map["O2+"]->nvib=1;
    name_map["O2+"]->theta[0]=2.74120e+03;

    // OH
    name_map["OH"]->cfs=2.5;
    name_map["OH"]->mol_weight=17.00800;
    name_map["OH"]->nvib=1;
    name_map["OH"]->theta[0]=5.37820e+03;

    // Si
    name_map["Si"]->cfs=1.5;
    name_map["Si"]->mol_weight=28.08550;
    name_map["Si"]->nvib=0;

    // SiO
    name_map["SiO"]->cfs=2.5;
    name_map["SiO"]->mol_weight=44.08550;
    name_map["SiO"]->nvib=1;
    name_map["SiO"]->theta[0]=1.78640e+03;

    // electron
    name_map["e"]->cfs=1.5;
    name_map["e"]->mol_weight=0.00055;
    name_map["e"]->nvib=0;

    for (ii=0; ii < SS.size(); ii++) {
        // check nvib was initialized for all species
        if (name_map[SS[ii]]->nvib == -1) {
            std::cout << name_map[SS[ii]]->nvib << std::endl;
            throw CanteraError("Error in StatMech.cpp",
                               "nvib not initialized!. \n\n");

        } else {
            // check that theta is initialized
            for (int i=0; i<name_map[SS[ii]]->nvib; i++) {
                if (name_map[SS[ii]]->theta[i] <= 0.0) {
                    throw CanteraError("Error in StatMech.cpp",
                                       "theta not initialized!. \n\n");
                }
            }

            // check that no non-zero theta exist
            // for any theta larger than nvib!
            for (int i=name_map[SS[ii]]->nvib; i<5; i++) {
                if (name_map[SS[ii]]->theta[i] != 0.0) {
                    std::string err = "bad theta value for "+SS[ii]+"\n";
                    throw CanteraError("StatMech.cpp",err);
                }
            } // done with for loop
        }

        // check mol weight was initialized for all species
        if (name_map[SS[ii]]->mol_weight == -1) {
            std::cout << name_map[SS[ii]]->mol_weight << std::endl;
            throw CanteraError("Error in StatMech.cpp",
                               "mol_weight not initialized!. \n\n");

        }

        // cfs was initialized for all species
        if (name_map[SS[ii]]->cfs == -1) {
            std::cout << name_map[SS[ii]]->cfs << std::endl;
            throw CanteraError("Error in StatMech.cpp",
                               "cfs not initialized!. \n\n");

        }

    } // done with sanity checks

    // mark it zero, dude
    return 0;
}

void StatMech::updateProperties(const doublereal* tt,
                                doublereal* cp_R, doublereal* h_RT,
                                doublereal* s_R) const
{

    std::map<std::string,species*>::iterator it;

    // get species name, to gather species properties
    species* s;

    // pointer to map location of particular species
    if (name_map.find(sp_name)  != name_map.end()) {
        s = name_map.find(sp_name)->second;
    } else {
        throw CanteraError("StatMech.cpp",
                           "species properties not found!. \n\n");
    }

    // translational + rotational specific heat
    doublereal ctr = 0.0;
    double theta = 0.0;

    // 5/2 * R for molecules, 3/2 * R for atoms
    ctr += GasConstant * s->cfs;

    // vibrational energy
    for (int i=0; i< s->nvib; i++) {
        theta = s->theta[i];
        ctr += GasConstant * theta * (theta* exp(theta/tt[0])/(tt[0]*tt[0]))/((exp(theta/tt[0])-1) * (exp(theta/tt[0])-1));
    }

    // Cp = Cv + R
    doublereal cpdivR = ctr/GasConstant + 1;

    // ACTUNG: fix enthalpy and entropy
    doublereal hdivRT = 0.0;
    doublereal sdivR  = 0.0;

    // return the computed properties in the location in the output
    // arrays for this species
    cp_R[m_index] = cpdivR;
    h_RT[m_index] = hdivRT;
    s_R [m_index] = sdivR;
}

void StatMech::updatePropertiesTemp(const doublereal temp,
                                    doublereal* cp_R, doublereal* h_RT,
                                    doublereal* s_R) const
{
    double tPoly[1];
    tPoly[0]  = temp;
    updateProperties(tPoly, cp_R, h_RT, s_R);
}

void StatMech::reportParameters(size_t& n, int& type,
                                doublereal& tlow, doublereal& thigh,
                                doublereal& pref,
                                doublereal* const coeffs) const
{
    species* s;

    n = m_index;
    type = STAT;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    for (int i = 0; i < 9; i++) {
        coeffs[i] = 0.0;
    }
    doublereal temp = coeffs[0];
    coeffs[1] = m_lowT;
    coeffs[2] = m_highT;

    // get species name, to gather species properties
    // pointer to map location of particular species
    if (name_map.find(sp_name)  != name_map.end()) {
        s = name_map.find(sp_name)->second;
    } else {
        throw CanteraError("StatMech.cpp",
                           "species properties not found!. \n\n");
    }

    double theta = 0.0;
    doublereal cvib = 0;

    // vibrational energy
    for (int i=0; i< s->nvib; i++) {
        theta = s->theta[i];
        cvib += GasConstant * theta * (theta* exp(theta/temp)/(temp*temp))/((exp(theta/temp)-1) * (exp(theta/temp)-1));
    }

    // load vibrational energy
    coeffs[3] = GasConstant * s->cfs;
    coeffs[4] = cvib;

}

void StatMech::modifyParameters(doublereal* coeffs)
{
}

}
