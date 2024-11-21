//
// Created by Liang on 2023/4/21.
//

#ifndef BCS_DATA_POST_ANALYSIS_H
#define BCS_DATA_POST_ANALYSIS_H
#include "Sim1D.h"
#include "StFlow.h"
#include "Boundary1D.h"
#include <stdio.h>
#include <fstream>
#include "cantera/base/Solution.h"
#include <filesystem>

namespace Cantera
{/*this file is designed to do sensitivity analysis for liquid-pool counterflow flame results*/
    class Data_post_analysis {
    public:
        Data_post_analysis(){
            std::cout<<"\n..............................................................................\n";
            std::cout<<"*************************Data post analysis is called*************************\n";
            std::cout<<"..............................................................................\n";}
        //For gas phase counterflow flame analysis
        void Initialize(Sim1D &flame, StFlow &flow, shared_ptr<Solution> sol, std::string sol_id);
        //For liquid pool conterflow flame analysis
        void Initialize(Sim1D &flame, StFlow &flow,Inlet1D_new &inlet, Inlet1D_new &outlet, shared_ptr<Solution> sol, std::string sol_id);

        //Restore flame with existed solution
        void Restore(const std::string &filename, const std::string &file_id )
        {m_flame->restore(filename, file_id);
            std::cout<<m_flame->nDomains()<<std::endl;}

        //print temperature for all points
        void Temperature(){for(size_t n=0; n<points;n++){std::cout<<solution_vector[n*n_variable+c_offset_T]<<std::endl;}}

        //return temeprature at point j
        double Temperature(size_t j){return solution_vector[j*n_variable+c_offset_T];}

        //i starts from the first specie, j is the loc_point
        double MassFraction(size_t i, size_t j){return solution_vector[i+j*n_variable+c_offset_Y];}

        //Output solution vector and molefraction to file
        void Output_Solution(std::string Result_file_name);

        //Output solution at liquid-gas phase
        void Output_liquid_pool(std::string Liquid_pool_solution_name);

        //Output reaction rate and constant to file
        void Output_Reaction_rate(std::string Reaction_rate_file_name);
        
        //Output heat_prod_rate to file
        void Output_Heat_Prod_Rate(std::string Heat_Prod_Rate_filename);
        
        //calculate sensitivity
        void get_temperature_reaction_sensitivity(size_t point_loc);
        void get_oxygen_sensitivity(size_t point_loc);
        void get_sensitivity(size_t point_loc, std::string variable);
        void get_sensitivity_all(std::string variable);

        void Output_Sensitivity(std::vector<std::string> var_vec, std::string Sensitivity_file_name);

        void perturb(size_t i, double dp ){m_flow->kinetics().setMultiplier(i,1+dp);}

        vector<double> solve_adjoint( size_t n_params, vector<double> dgdx, double(*g_ptr)()=nullptr, double dp=1e-5 );
        
        void solve_adjoint_all(size_t n_params, Array2D* dgdx_all, double(*g_ptr)()=nullptr, double dp=1e-5 );

        void check_path(std::string path_name){
            std::filesystem::path folderPath(path_name);
            if (std::filesystem::exists(folderPath)) {
            } else {
                std::filesystem::create_directories(folderPath);
            }
        }
    protected:
        Sim1D *m_flame;
        StFlow *m_flow;
        Inlet1D_new *m_inlet, *m_outlet;
        Phase_liquid *m_pool;
        shared_ptr<Solution> m_sol;
        vector<double> r_Mean_mw,solution_vector,Sens_reaction;
        std::vector<std::string> reaction_name_vector;
        size_t nsp, points, n_variable, n_reaction;
        Array2D dfdp, Sens_reaction_all;
        bool dfdp_flag;
        std::string folder_name;
    };
}
#endif //BCS_DATA_POST_ANALYSIS_H
