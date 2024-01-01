//
// Created by Liang on 2023/3/27.
//

#include <cantera/base/stringUtils.h>
#include "cantera/oneD/Phase_liquid.h"
#include <stdio.h>

namespace Cantera
{
    //constructor
    Phase_liquid::Phase_liquid(): sum_mole_weight_inlet(0), sum_mole_weight_outlet(0), temperature_inlet(298.15), temperature_outlet(298.15), nsp_liquid(0)
    {   //mole_fraction_inlet(100,0),mole_fraction_outlet(100,0),liquid_enthalpy(100,0),delta_H_sp(100)
        mole_fraction_inlet.resize(name.size(),0);
        mole_fraction_outlet.resize(name.size(),0);
        mole_fraction_gas.resize(name.size(),0);
        liquid_enthalpy.resize(name.size(),0);
        delta_H_sp.resize(name.size(),0);
        for (size_t n=0; n<name.size();n++){m_speciesIndices[name[n]]=n;m_speciesName[n]=name[n];}
    }

    void Phase_liquid::setMoleFraction(std::string name, const double mole)//initialization phase by molarfraction
    {
            size_t m=speciesIndex(name);
        if (abs(mole)>1e-10){
            Loc_sp.push_back(m);
            mole_fraction_inlet[m]=mole;
            mole_fraction_outlet[m]=mole;
            nsp_liquid++;}
            update();
    }

    void Phase_liquid::setMoleFraction(const std::string& fuel_comp)//initialization phase by molarfraction
    {
        for(const auto& sp: parseCompString(fuel_comp)){
            size_t loc=speciesIndex(sp.first);
            if (loc==npos){throw CanteraError("Phase_Liquid::setMoleFraction", "Unknown species'{}'", sp.first);}
            if (abs(sp.second)>1e-10){
            Loc_sp.push_back(loc);
            mole_fraction_inlet[loc]=sp.second;
            mole_fraction_outlet[loc]=sp.second;
            nsp_liquid++;}
        }
        update();
    }

    void Phase_liquid::Build_interface(double *xb, StFlow *m_flow ){
        //copy molar fractions of corresponding activated species in gas phase
        for (size_t i = 0; i < nsp_liquid; i++) {mole_fraction_gas[Loc_sp[i]] = m_flow->molar_fraction(xb, m_flow->componentIndex(name[Loc_sp[i]]) - 5, 0);}

        if (nsp_liquid==1) {
            size_t loc=Loc_sp[0];
            temperature_interface_gas= SatTemperature(name[loc],m_flow->pressure()*mole_fraction_gas[loc]/mole_fraction_outlet[loc]);}

        else if (nsp_liquid==2) {
            vector<double> T = {0, 300};//T0: last obtained temperature, T1: new obtained temperature
            size_t loc0=Loc_sp[0];
            size_t loc1=Loc_sp[1];

            T[1] = SatTemperature(name[loc1], m_flow->pressure() * mole_fraction_gas[loc1] / mole_fraction_inlet[loc1]);
            while (abs(T[1] - T[0]) > 1e-10 || T[0] == 0) {
                double mole_ratio=mole_fraction_gas[loc0] / mole_fraction_gas[loc1];
                double Psat_ratio=SatPressure(name[loc1], T[1]) / SatPressure(name[loc0], T[1]);
                mole_fraction_outlet[loc1] = 1.0 / (1.0 + mole_ratio*Psat_ratio);
                mole_fraction_outlet[loc0] = 1.0 - mole_fraction_outlet[loc1];

                T[0] = T[1];
                T[1] = SatTemperature(name[loc1], m_flow->pressure() * mole_fraction_gas[loc1] / mole_fraction_outlet[loc1]);}
            temperature_interface_gas = T[1];
        }
        else {throw CanteraError("Phase_Liquid::Build_interface", "Number of species cannot be larger than 2");}
        update();
    }


    void Phase_liquid::update(){//update sum_molecular weight and sum_liquid enthalpy and temperature
        sum_mole_weight_inlet=0;
        sum_mole_weight_outlet=0;
        sum_liquid_ethalpy=0;
        //std::cout<<"before sum_mole_weight="<<sum_mole_weight<<std::endl;
        for(size_t n=0; n < nsp_liquid; n++){//update for temperature change
            liquid_enthalpy[Loc_sp[n]] = H_formation[Loc_sp[n]] + Cp[Loc_sp[n]] * (temperature_outlet-298.15);}
        //update for adding species
        for(size_t n=0;n < nsp_liquid; n++){//std::cout<<n<<" sum_mole_weight="<<sum_mole_weight;
            sum_mole_weight_inlet += mole_fraction_inlet[Loc_sp[n]] * molecular_weight[Loc_sp[n]];
            sum_mole_weight_outlet += mole_fraction_outlet[Loc_sp[n]] * molecular_weight[Loc_sp[n]];
            sum_liquid_ethalpy += mole_fraction_inlet[Loc_sp[n]] * liquid_enthalpy[Loc_sp[n]];}
        //std::cout<<std::endl<<"after sum_mole_weight="<<sum_mole_weight<<mole_fraction[0]*molecular_weight[0]<<std::endl;
    }

    //![J/kmol] liquid enthalpy of each single component //C7H16 J/kmol @298K/https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=2
    double Phase_liquid::Liquid_enthpy(std::string name, double temperature){
        setTemperature_outlet(temperature);
        return liquid_enthalpy[speciesIndex(name)];//[J/kmol]
    }

    //![J/kmol] enthalpy of component in gas phase
    double Phase_liquid::enthpy_gas(std::string name, StFlow *m_flow){
        return m_flow->Enthalpy()[m_flow->componentIndex(name)-c_offset_Y];
    }

    //![J/kmol] enthalpy change between liquid and gas in molar basis
    double Phase_liquid::Delta_H_mole(doublereal *xb, StFlow *m_flow){
        setTemperature_outlet(xb[c_offset_T]);
        Build_interface(xb, m_flow);
        delta_H_sum=0;
        for (size_t n=0;n < nsp_liquid; n++){
        delta_H_sum += mole_fraction_outlet[Loc_sp[n]] * Vap_enthpy(name[Loc_sp[n]],temperature_outlet);//use vaporization enthalpy directly
        //std::cout<<name[Loc_sp[n]]<<" "<<mole_fraction_inlet[Loc_sp[n]]<<" "<<mole_fraction_outlet[Loc_sp[n]]<<" "<<mole_fraction_gas[Loc_sp[n]]<<" "<<temperature_interface_gas<<std::endl;
        }
        return delta_H_sum;
    }
    
    double Phase_liquid::Delta_H_mole(double t){
        setTemperature_outlet(t);
        delta_H_sum=0;
        for (size_t n=0;n < nsp_liquid; n++){
        delta_H_sum+=mole_fraction_outlet[Loc_sp[n]]*Vap_enthpy(name[Loc_sp[n]],temperature_outlet);//use vaporization enthalpy directly
        }
        return delta_H_sum;
    }

    //![K] return temerpature at the liquid-gas interface
    double Phase_liquid::Temperature_interface(double *xb, StFlow *m_flow){ Build_interface(xb, m_flow);
        return temperature_interface_gas;}

    //*[K] calculate temperature based on average Tsat of each specie
    /*double Phase_liquid::Temperature_interface(double *xb, StFlow *m_flow){
        temperature_interface_gas=0;
        for (size_t i = 0; i < name.size(); i++) {mole_fraction_gas[i] = m_flow->molar_fraction(xb, m_flow->componentIndex(name[i]) - 5, 0);}
        for(size_t i=0;i<nsp_liquid;i++){
        temperature_interface_gas+= 0.5*SatTemperature(name[i],m_flow->pressure()*mole_fraction_gas[i]/mole_fraction_outlet[i]);
        }
        //std::cout<<"Temperature"<<temperature_interface_gas<<std::endl;
         return temperature_interface_gas;}*/
         
//https://app.knovel.com/kn/resources/kt002UT9T2/kpYHTPPCC4/eptble/itable?b-toc-cid=kpYHTPPCC4&b-toc-title=Yaws%27+Handbook+of+Thermodynamic+and+Physical+Properties+of+Chemical+Compounds&b-toc-url-slug=enthalpy-vaporization&columns=6%2C1%2C3%2C4%2C5%2C7%2C8%2C9%2C10%2C11%2C13%2C14&q=heptane

    double Phase_liquid::Vap_enthpy(std::string name, double temperature){
        double C, n, Tcr;
        if(name=="NC7H16" || name=="C7H16"){
            C=49.73;
            n=0.386;
            Tcr=540.26;}
        else if(name=="C2H5OH"){
            C=43.122;
            n=0.079;
            Tcr=516.25;}
        else{throw CanteraError("Phase_liquid::Vap_enthalpy","unknown species '{}' ",name);}

        double delta_H_sum=C*std::pow((1-temperature/Tcr),n)*1e6;
        return delta_H_sum;
    }

    double Phase_liquid::SatPressure(std::string name, double temperature) {//http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Ethanol
        double Tc = temperature - 273.15; //converted from Kelvin to C
        doublereal c_1, c_2, c_3, c_4, c_5;
        double P_Hg, P_pa;
        if (name == "C2H5OH") {
            if (temperature>516 || temperature<159){ std::cout << temperature << "K is out of Ethanol range" << std::endl;abort(); }
            c_1 = 23.844;
            c_2 = -2864.2;
            c_3 = -5.0474;
            c_4=3.7448E-11;
            c_5=2.7361E-07;
        }
        else if (name=="NC7H16" || name=="C7H16") {
            if (temperature > 540 || temperature < 182) {std::cout << temperature << "K is out of Heptane range" << std::endl; }
            c_1 = 65.026;
            c_2 = -3818.8;
            c_3 = -21.684;
            c_4=1.0387E-02;
            c_5=1.0206E-14;
        }
        else{throw CanteraError("Phase_liquid::SatPressure","unknown species '{}' ",name);}
        P_Hg = std::pow(10, c_1+c_2/temperature+c_3*std::log10(temperature)+c_4*temperature+c_5*temperature*temperature);//pressure in mmHg
        P_pa = P_Hg / 760 * 101325;//conversion from mmHg to Pa
        return P_pa;
    }

    double Phase_liquid::SatTemperature(std::string name, double pressure){
        double T1=273, T2=550, T_mid;//T2 depends on upper limit of valid temperature
        while(std::abs(T1-T2)>1e-12){
            T_mid=0.5*(T1+T2);
            pressure>SatPressure(name, T_mid)? T1=T_mid:T2=T_mid;}

        if (std::abs(T1-273.15)<0.1 || std::abs(T1-550)<0.1){throw CanteraError("Phase_liquid::SatTemperature","'{}' K is out of '{}' range:  ", T1, name);}

        return T1;
    }
};

