/**
 *  @file AVBPTransport.cpp
 *  Simplified AVBP transport properties for ideal gas mixtures.
 */

/* $Author: B. Franzelli (v. 1.7) $
 * $Revision: A. Felden (v 2.1-2.3) $
 * $Date: 01/2018 $
 */

#include "cantera/transport/AVBPTransport.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/Parser.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

namespace Cantera
{

// AVBPTransport::AVBPTransport() :
//     m_lambda(0.0)
// {
// }

void AVBPTransport::init(ThermoPhase* thermo, int mode, int log_level)
{

    GasTransport::init(thermo, mode, log_level);

    ifstream avbp_mixturedb("mixture_database.dat");

    if (avbp_mixturedb.good()){                                                                                 
    	  //cout << thermo->name() << thermo->id() << endl;
    	  //std::string  mixture_id; 
          //mixture_id = thermo->name();
          //cout<<"INFO: using simplified transport data for "<< m_thermo->name()<<" from the 'mixture_database.dat' file."<<endl;

	  AVBPTransport::read_mixture("mixture_database.dat");

    }
    else {                                                                                 
          cout<<"FATAL ERROR: cannot find required input files (mixture_database.dat) for simplified transport model !"<<endl;                                  
	  exit(-1);
    }                                                                            



//    // Print checks
//    cout<<endl<<"++++++++++AVBP Transports is used:++++++++++++"<<endl;
//    cout<<"-Reference viscosity: "<<avbp_mu0<<endl;
//    cout<<"-Reference temperature: "<<avbp_T0<<endl;
//    if(avbp_beta<0){
//      cout<<"-Exponent for power law: "<<-avbp_beta<<endl;
//    }
//    else if(avbp_beta>0){
//       cout<<"-Exponent for Second Sutherland Constant: "<<avbp_beta<<endl;
//    }
//    else{
//       cout<<"The molecular viscosity used is temperature indipendent and equal to: "<<avbp_mu0;
//    }
//    cout<<"-Prandtl number: "<<avbp_Prandtl<<endl;
//    cout<<"-Schmidt and Lewis number for: "<<endl;
//    for(int i=0;i<m_nsp;i++){
//      cout<<"*Specie "<<i<<": "<<avbp_Sch[i]<<" "<<avbp_Le[i]<<endl;
//    }

    double m_cp = 0.0;
    double m_density = 0.0;
    m_cp = m_thermo->cp_mass();
    m_density = m_thermo->density();
    if(avbp_beta<0){
      double avbp_absbeta;
      avbp_absbeta = - avbp_beta;            
      m_viscmix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }
    else if(avbp_beta>0){
      double coeff;
      coeff= (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5); 
      m_viscmix = avbp_mu0 * coeff * pow(m_temp, 1.5) / (m_temp + avbp_beta);  
    }
    else{
    m_viscmix = avbp_mu0;
    }
    m_lambda = m_viscmix * m_cp / avbp_Prandtl;
}


double AVBPTransport::viscosity()
{
    update_T();
    update_C();

    double vismix = 0.0;

    // AVBP Transport Properties
    if(avbp_beta<0){
          double avbp_absbeta;
          avbp_absbeta = - avbp_beta;
          vismix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }
    else if(avbp_beta>0){
          double coeff;
          coeff= (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5);
          vismix = avbp_mu0 * coeff * pow(m_temp, 1.5) / (m_temp + avbp_beta);
    }
    else{
        vismix = avbp_mu0;
    }   
     
    m_viscmix = vismix;
    return m_viscmix;
}


double AVBPTransport::thermalConductivity()
{
    update_T();
    update_C();

    double m_cp = 0.0;
    double m_vismix = 0.0;
    m_cp = m_thermo->cp_mass();
    
    if(avbp_beta<0){
          double avbp_absbeta;
          avbp_absbeta = - avbp_beta; 
          m_vismix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }else if(avbp_beta>0){   
          double coeff;
          coeff = (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5); 
          m_vismix = avbp_mu0 * coeff * pow (m_temp, 1.5) / (m_temp + avbp_beta);
    }else{
          m_vismix = avbp_mu0;
    }

    m_lambda = m_vismix * m_cp / avbp_Prandtl;
    return m_lambda;
}


void AVBPTransport::getMixDiffCoeffs(double* const d)
{
    update_T();
    update_C();

    double m_cp = 0.0;
    double m_density = 0.0;
    m_cp = m_thermo->cp_mass();
    m_density = m_thermo->density();
 
    double m_vismix = 0.0;
    if(avbp_beta<0){
          double avbp_absbeta;
          avbp_absbeta = - avbp_beta;
          m_vismix = avbp_mu0 * pow(m_temp/avbp_T0,avbp_absbeta);
    }else if(avbp_beta>0){      
          double coeff;
          coeff = (avbp_T0 + avbp_beta) / pow(avbp_T0,1.5);
          m_vismix = avbp_mu0 * coeff * pow(m_temp,1.5)/ (m_temp + avbp_beta);
    }else{
          m_vismix = avbp_mu0;
    }
 
    for (size_t k=0; k<m_nsp; k++){
          d[k]=m_vismix / m_density / avbp_Sch[k];
    }
}

/**
 *  @internal This is called whenever a transport property is
 *  requested from ThermoSubstance if the temperature has changed
 *  since the last call to update_T.
 */
void AVBPTransport::update_T()
{
    double t = m_thermo->temperature();
    if (t == m_temp && m_nsp == m_thermo->nSpecies()) {
        return;
    }
    if (t <= 0.0) {
        throw CanteraError("AVBPTransport::update_T",
                           "negative temperature {}", t);
    }

    GasTransport::update_T();

}

void AVBPTransport::update_C()
{
    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_thermo->getMoleFractions(m_molefracs.data());

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}


/**
 *  Parse the specific mixture_param.dat input file of 
 *  AVBP v7 ++. Use the Parser.cpp located in ../base
 */
void AVBPTransport::read_mixture(std::string inputfile)
{
    Parser parser;
    Param* param;

    parser.parseFile(inputfile);
    // Get number of mixtures in database
    size_t n_mixtures = parser.nbParamOccurence("mixture_name"); 
    size_t idx_beg = std::string::npos;
    size_t idx_end = std::string::npos;
    // Loop through all occurences of keyword "mixture name"
    for (size_t i=0; i< n_mixtures; ++i) {
      param = parser.getParam("mixture_name",i);
      std::string current_str = param->readString(0);
      // Check if mixture_name is the one requested in the cti file
      if ( current_str == m_thermo->name() ) {
        // Store bounding idx to get important info
	idx_beg = parser.getParamNumber("mixture_name",i);
	idx_end = parser.getParamNumber("mixture_name",i+1);
      }
    }  
    // Get an error if requested mixture does not match any entry in database
    if (idx_beg == idx_end) {
    	        cout<<"FATAL ERROR: cannot find the requested mixture in the 'mixture_database.dat' file.\n \
    	        Make sure that the name of your gas instance in the Cantera mechanism (.cti/.xml) file matches that provided in the 'mixture_database.dat'."<<endl;                                  
    	        exit(-1);
    }
    // Otherwise read important data. NB if keyword does not exist, error is managed in Parser.cpp 
    std::string avbp_law ;
    // PRANDTL
    param = parser.getParam("prandtl_number", idx_beg, idx_end);
    avbp_Prandtl = param->readDouble(0);
    // SCH LE
    param = parser.getParam("species_Schmidt_number", idx_beg, idx_end);
    avbp_Le.resize(m_nsp);
    avbp_Sch.resize(m_nsp);
    for (size_t i=0; i<m_nsp; ++i) {
	    avbp_Sch[i] = param->readDouble(i); 
      avbp_Le[i]= avbp_Sch[i] / avbp_Prandtl;
    }
    // VISCO
    param = parser.getParam("mu_ref", idx_beg, idx_end);
    avbp_mu0 = param->readDouble(0);
    param = parser.getParam("T_ref", idx_beg, idx_end);
    avbp_T0 = param->readDouble(0);
    param = parser.getParam("viscosity_law_coeff", idx_beg, idx_end);
    avbp_beta = param->readDouble(0);
    param = parser.getParam("viscosity_law", idx_beg, idx_end);
    avbp_law = param->readString(0);
    if (avbp_law == "power") {
	    avbp_beta = - avbp_beta;
    }
}

double AVBPTransport::pressure_ig()
{
    return (m_thermo->molarDensity() * GasConstant *
            m_thermo->temperature());
    
}

}