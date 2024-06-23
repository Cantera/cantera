//
// Created by Liang on 2023/1/8.
//

#ifndef BCS_PHASE_LIQUID_H
#define BCS_PHASE_LIQUID_H
#include "cantera/thermo/Phase.h"
#include "StFlow.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
namespace Cantera
{/**The class of boundary for liquid pool used in one-dimensional spatial domains.
 * setup inlet molar fraction of liquid pool and get outlet conditions
*/
class Phase_liquid
{
public:
    //!constructor
    Phase_liquid();
    
    //![K] set inlet temperature for liquid_phase, default temperature is 298.15K
    void setTemperature_inlet(double t){temperature_inlet=t;};
    
    //!return inlet temeprature for liquid_pool
    double Temperature_inlet(){return temperature_inlet;};

    //![K] set outlet temperature for liquid_phase, default temperature is 298.15K
    void setTemperature_outlet(double t){temperature_outlet=t;
       update();};

    //![K] return inlet temperature for liquid_phase
    double Temperature_outlet(){return temperature_outlet;};

    //!set inlet mole fraction of liquid pool
    void setMoleFraction(std::string name, const double mole);

    //!set inlet mole fraction of liquid pool
    void setMoleFraction(const std::string& fuel_composition);

    //!calculate conditions (temperature, molar fraction) at liquid side of interface based on gas phase (support for mixture less than two components)
    void Build_interface(double *xb, Cantera::StFlow *m_flow );

    //!return inlet mole fraction by number
    double MoleFraction_inlet(int m){return mole_fraction_inlet[m];}
    //!return inlet mole fraction by name
    double MoleFraction_inlet(std::string species_name){return MoleFraction_inlet(speciesIndex(species_name));}

    //!return outlet mass fraction by number
    double MassFraction_inlet(int m) {return mole_fraction_inlet[m]*molecular_weight[m]/sum_mole_weight_inlet;}
    //!return outlet mass fraction by name
    double MassFraction_inlet(std::string species_name){return MassFraction_inlet(speciesIndex(species_name));}

    //!return outlet mole fraction by number
    double MoleFraction_outlet(int m){return mole_fraction_outlet[m];}
    //!return outlet mole fraction by name
    double MoleFraction_outlet(std::string species_name){ return MoleFraction_outlet(speciesIndex(species_name));}

    //!return outlet mass fraction by number
    double MassFraction_outlet(int m){return mole_fraction_outlet[m]*molecular_weight[m]/sum_mole_weight_outlet;};
    //!return outlet mass fraction by name
    double MassFraction_outlet(std::string species_name){ return MassFraction_outlet(speciesIndex(species_name));};

    //!update sum_molecular weight and sum_liquid enthalpy and temperature
    void update();

    //![J/kmol] liquid enthalpy of each single component at specified temperature
    double Liquid_enthpy(std::string name, double temperature);

    //![J/kmol] liquid enthalpy of mixture
    double Liquid_enthpy(){return sum_liquid_ethalpy;}

    //![J/kmol] enthalpy of component in gas phase
    double enthpy_gas(std::string name, StFlow *m_flow);

    //![J/kmol] enthalpy change between liquid and gas in molar basis
    double Delta_H_mole(doublereal *xb, StFlow *m_flow);
    double Delta_H_mole(double t);

    //![J/kg] enthalpy change at outlet from liquid to gas in mass basis
    double Delta_H_mass(doublereal *xb,StFlow *m_flow){return Delta_H_mole(xb,m_flow)/sum_mole_weight_outlet;}
    double Delta_H_mass(double t){return Delta_H_mole(t)/sum_mole_weight_outlet;}

    //![K] return temerpature at the liquid-gas interface
    double Temperature_interface(double *xb, StFlow *m_flow);

    //![J/kmol] vaporization enthalpy of each single component at specified temperature
    double Vap_enthpy(std::string name, double temperature);

    //![Pa] return of saturated pressure of specie (name) at temperature [K]
    double SatPressure(std::string name, double temperature);

    //![K] return of saturated temperature of specie (name) at pressure [Pa]
    double SatTemperature(std::string name, double pressure);

    //!return the number of species in liquid pool
    size_t Nsp_liquid(){return nsp_liquid;}

    //!return location of species at liquid phase
    size_t speciesIndex(const std::string nameStr){
        size_t loc=npos;
        auto it=m_speciesIndices.find(nameStr);
        if (it!=m_speciesIndices.end()){return it->second;}
        return loc;
    }
    //!return name of species at liquid phase
    std::string speciesName(const size_t n){return m_speciesName.find(n)->second;}

protected:
    vector<double> mole_fraction_inlet, mole_fraction_outlet, mole_fraction_gas;
    vector<double> liquid_enthalpy,delta_H_sp;//[J/kmol]
    std::vector<size_t> Loc_sp;//location of activated species in pool

    std::vector<std::string> name{"C2H5OH", "NC7H16", "C7H16", "H2O"};
    vector<double>    H_formation{-276*1e6, -224*1e6, -224*1e6, -285.83*1e6};//J/kmol
    vector<double>             Cp{112*1e3, 224.64*1e3, 224.64*1e3, 75.37*1e3};//J/kmol/K}
    vector<double> molecular_weight{46.07, 100.21, 100.21, 18.015};//kg/kmol

    //!map of species names to indices
    std::map<std::string, size_t> m_speciesIndices;
    //!map of species indices to names
    std::map<size_t, std::string> m_speciesName;
    
    double sum_mole_weight_inlet, sum_mole_weight_outlet, sum_liquid_ethalpy,temperature_inlet, temperature_outlet, temperature_interface_gas;
    double delta_H_sum;

    size_t nsp_liquid ;//number of liquid species activated
};
}
#endif //BCS_PHASE_LIQUID_H
