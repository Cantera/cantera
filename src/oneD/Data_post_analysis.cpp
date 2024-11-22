//
// Created by Liang on 2023/4/21.
//

#include "cantera/oneD/Data_post_analysis.h"

namespace Cantera{

    void Data_post_analysis::Initialize(Sim1D &flame, StFlow &flow, shared_ptr<Solution> sol, std::string sol_id){
        auto gas(sol->thermo());
        m_flame=&flame;
        m_flow=&flow;
        m_sol=sol;
        dfdp_flag=false;
        //assgined the number of points to points
        points=m_flame->domain(1).nPoints();//std::cout<<"Points:"<<points<<"\n";
        //assgined the number of components to nsp (not include first 5 variables)
        nsp=gas->nSpecies();
        //the number of variables
        n_variable=nsp+c_offset_Y;
        //get solution vector from flame
        solution_vector.resize(n_variable*points);
        for(size_t n=0;n<flame.size();n++){
            solution_vector[n]=(flame.value(1,n%n_variable,n/n_variable));
        }
        //the number of reactions in kinetics
        n_reaction=m_flow->kinetics().nReactions();
        //get reaction name vector from flow
        reaction_name_vector.resize(n_reaction);
        for(size_t n=0;n<n_reaction;n++){reaction_name_vector[n]=m_flow->kinetics().reactionString(n);}

        folder_name="Post_process/"+sol_id+"K";
        check_path(folder_name);
        Output_Solution(folder_name+"/Results.txt");
        Output_Reaction_rate(folder_name+"/Reaction_rate.txt");
        Output_Heat_Prod_Rate(folder_name+"/Heat_prod_rate.txt");
        Output_Sensitivity({"T"},folder_name+"/Sensitivity.txt");
        std::cout<<"Data post analysis has been finished\n";
    
    }
    
    void Data_post_analysis::Initialize(Sim1D &flame, StFlow &flow,Inlet1D_new &inlet, Inlet1D_new &outlet, shared_ptr<Solution> sol, std::string sol_id){
        Initialize(flame, flow, sol, sol_id);
        m_inlet=&inlet;
        m_outlet=&outlet;
        m_pool=m_inlet->Pool();
        Output_liquid_pool(folder_name+"/Result_pool.txt");
    }

    void Data_post_analysis::Output_Solution(std::string Result_file_name){
        std::ofstream Result(Result_file_name);
        //print result to Result.txt
        for(size_t m=0;m<points;m++){
            //print title to Result file
            if(m==0){
                for(size_t n=0;n<nsp+n_variable;n++){
                    //print grid to first column
                    if(n==0){Result<<std::setw(18)<<std::left<<"Grid[m]";}
                    std::string component_name=m_flame->domain(1).componentName(n);
                    //print title of molar fraction
                    if(n>=c_offset_Y&&n<n_variable){component_name="Y_"+m_flame->domain(1).componentName(n);}
                    //print title of mass fraction
                    if(n>=n_variable) {component_name="X_"+m_flame->domain(1).componentName(n-n_variable+c_offset_Y);}
                    Result<<std::setw(18)<<std::left<<component_name;}
                Result<<std::endl;}

            //print variables values to Result file
            r_Mean_mw.push_back(0);
            for(size_t n=0;n<nsp+n_variable;n++)
            {// 1st column: print grid value to first column
                if(n==0){Result<<std::setw(18)<<std::left<<std::scientific<<m_flame->domain(1).grid(m);}
                //From the 2nd column to n_variable column: mass fraction range
                if(n<n_variable){
                    //calculate the molecular weight at point m
                    if(n>=c_offset_Y){r_Mean_mw[m]+=MassFraction(n-c_offset_Y,m)/m_flow->Mw(n-c_offset_Y);}
                    //print molar fraction value
                    Result<<std::setw(18)<<std::left<<std::scientific<<m_flame->value(1,n,m);}
                //From n_variable column to last column: print molar fraction value
                if(n>=n_variable){Result<<std::setw(18)<<std::left<<std::scientific<<MassFraction(n-n_variable,m)/r_Mean_mw[m]/m_flow->Mw(n-n_variable);}
            }
            Result<<std::endl;
        }
    }

    void Data_post_analysis::Output_liquid_pool(std::string Liquid_pool_solution_name){
        std::ofstream Liquid_pool_solution(Liquid_pool_solution_name);
        std::vector<std::string> liquid_temp_title_vector{"Temperature_inlet","Temperature_outlet"};
        vector<double> liquid_temp_vector={m_pool->Temperature_inlet(),m_pool->Temperature_outlet()};

        for(size_t m=0;m<2;m++){//inlet and outlet
            vector<double> liquid_result_vector;
            Liquid_pool_solution<<std::setw(20)<<std::left<<liquid_temp_title_vector[m];
            size_t Num_liquid=m_pool->Nsp_liquid();
            for(size_t n=0;n<2*Num_liquid;n++){
                if(n<Num_liquid){
                    Liquid_pool_solution<<std::setw(20)<<std::left<<"X_"+m_pool->speciesName(n%Num_liquid);
                    liquid_result_vector.push_back(m==0? m_pool->MoleFraction_inlet(n%Num_liquid):m_pool->MoleFraction_outlet(n%Num_liquid));}
                else{Liquid_pool_solution<<std::setw(20)<<std::left<<"Y_"+m_pool->speciesName(n%Num_liquid);
                    liquid_result_vector.push_back(m==0? m_pool->MassFraction_inlet(n%Num_liquid):m_pool->MassFraction_outlet(n%Num_liquid));}
            }
            Liquid_pool_solution<<std::endl;
            for(size_t n=0;n<2*Num_liquid;n++){
                if(n==0){Liquid_pool_solution<<std::setw(20)<<std::left<<liquid_temp_vector[m];}
                Liquid_pool_solution<<std::setw(20)<<std::left<<liquid_result_vector[n];}
            Liquid_pool_solution<<std::endl;
        }
    }

    void Data_post_analysis::Output_Reaction_rate(std::string Reaction_rate_file_name){

        std::ofstream Reaction_Rate(Reaction_rate_file_name);
        Array2D forward_rate(n_reaction,points),reverse_rate(n_reaction,points),net_reaction_rate(n_reaction,points),forward_rate_coeff(n_reaction,points),reverse_rate_coeff(n_reaction,points),eq_rate_coeff(n_reaction,points);
        for(size_t j=0;j<points;j++){
            m_flow->setGas(&solution_vector[0],j);
            m_flow->kinetics().getFwdRatesOfProgress(&forward_rate(0,j));
            m_flow->kinetics().getRevRatesOfProgress(&reverse_rate(0,j));
            m_flow->kinetics().getNetRatesOfProgress(&net_reaction_rate(0,j));
            m_flow->kinetics().getFwdRateConstants(&forward_rate_coeff(0,j));
            m_flow->kinetics().getRevRateConstants(&reverse_rate_coeff(0,j));
            m_flow->kinetics().getEquilibriumConstants(&eq_rate_coeff(0,j));
        }
        //print reaction rate to Reaction_rate.txt
        for(size_t m=0;m<points;m++){
            //print title to Reaction_rate file
            if(m==0){
                for(size_t n=0;n<6*n_reaction;n++){
                    //print grid to first column
                    if(n==0){Reaction_Rate<<std::setw(50)<<std::left<<"Grid[m]";}
                    std::string reaction_name=reaction_name_vector[n%n_reaction];
                    if(n<n_reaction){reaction_name+="[f]";}
                    if(n>=n_reaction&&n<2*n_reaction){reaction_name+="[r]";}
                    if(n>=2*n_reaction&&n<3*n_reaction){reaction_name+="[net]";}
                    if(n>=3*n_reaction&&n<4*n_reaction){reaction_name+="[const_f]";}
                    if(n>=4*n_reaction&&n<5*n_reaction){reaction_name+="[const_r]";}
                    if(n>=5*n_reaction&&n<6*n_reaction){reaction_name+="[const_eq]";}
                    Reaction_Rate<<std::setw(50)<<std::left<<reaction_name;}
                Reaction_Rate<<std::endl;}
            //print reaction rate to Reaction_rate file
            for(size_t n=0;n<6*n_reaction;n++){
                // print grid value to first column
                if(n==0){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<m_flame->domain(1).grid(m);}
                //forward reaction rate
                if(n<n_reaction){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<forward_rate(n,m);}
                //reverse reaction rate
                if(n>=n_reaction&&n<2*n_reaction){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<reverse_rate(n-n_reaction,m);}
                //net reaction rate
                if(n>=2*n_reaction&&n<3*n_reaction){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<net_reaction_rate(n-2*n_reaction,m);}
                //forward reaction rate constant
                if(n>=3*n_reaction&&n<4*n_reaction){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<forward_rate_coeff(n-3*n_reaction,m);}
                //reverse reaction rate constant
                if(n>=4*n_reaction&&n<5*n_reaction){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<reverse_rate_coeff(n-4*n_reaction,m);}
                //equilibrium constant
                if(n>=5*n_reaction&&n<6*n_reaction){Reaction_Rate<<std::setw(50)<<std::left<<std::scientific<<eq_rate_coeff(n-5*n_reaction,m);}
            }
            Reaction_Rate<<std::endl;
        }
    }

    void Data_post_analysis::Output_Sensitivity(std::vector<std::string> var_vec, std::string Sensitivity_file_name){
    
        //creat vector of filename output
        size_t dot_pos=Sensitivity_file_name.find_last_of(".");
        std::vector<std::string> file_name_vec;
        std::string base_name=Sensitivity_file_name.substr(0,dot_pos);
        std::string ext_name=Sensitivity_file_name.substr(dot_pos);
        for(size_t n=0;n<var_vec.size();n++){
        file_name_vec.push_back(base_name+"_"+var_vec[n]+ext_name);
        std::cout<<"Sensitivity output of "<<var_vec[n]<<" is in "<<file_name_vec[n]<<"\n";}
         
        //print Sens_reaction to Sensitivity.txt
        /*for(size_t m=0;m<points;m++){
          for(size_t n=0;n<var_vec.size();n++){
            std::cout<<"Calculating sensitivity at point: "<<m<<"/"<<points<<"\n";
            get_sensitivity(m,var_vec[n]);
            //get_temperature_reaction_sensitivity(m);
            //get_oxygen_reaction_sensitivity(m);
            std::ofstream Sensitivity_file(file_name_vec[n],std::ios::app);
            if(m==0){//print title to Sensitivity file
                for(size_t n=0;n<n_reaction;n++){//print grid to first column
                    if(n==0){Sensitivity_file<<std::setw(50)<<std::left<<"Grid[m]"<<"\t";}
                    std::string reaction_name=reaction_name_vector[n%n_reaction];
                    Sensitivity_file<<std::setw(50)<<std::left<<reaction_name<<"\t";}
                Sensitivity_file<<std::endl;}
            //print sensitivity to Sensitivity file
            for(size_t n=0;n<n_reaction;n++){
                if(n==0){Sensitivity_file<<std::setw(56)<<std::left<<m_flame->domain(1).grid(m);}
                Sensitivity_file<<std::setw(56)<<std::left<<std::scientific<<Sens_reaction[n];}
            Sensitivity_file<<std::endl;
        }}*/
        for(size_t n=0;n<var_vec.size();n++){
            get_sensitivity_all(var_vec[n]);
            std::ofstream Sensitivity_file(file_name_vec[n],std::ios::app);
        for(size_t m=0;m<points;m++){
                std::cout<<"Outputing sensitivity of "+var_vec[n]+" at point: "<<m+1<<"/"<<points<<"\n";
                if(m==0){//print title to Sensitivity file
                    for(size_t n=0;n<n_reaction;n++){//print grid to first column
                        if(n==0){Sensitivity_file<<std::setw(50)<<std::left<<"Grid[m]"<<"\t";}
                        std::string reaction_name=reaction_name_vector[n%n_reaction];
                        Sensitivity_file<<std::setw(50)<<std::left<<reaction_name<<"\t";}
                    Sensitivity_file<<std::endl;}
                //print sensitivity to Sensitivity file
                for(size_t n=0;n<n_reaction;n++){
                    if(n==0){Sensitivity_file<<std::setw(56)<<std::left<<m_flame->domain(1).grid(m);}
                    Sensitivity_file<<std::setw(56)<<std::left<<std::scientific<<Sens_reaction_all.value(n,m);}
                Sensitivity_file<<std::endl;
            }}
    }
    
    void Data_post_analysis::Output_Heat_Prod_Rate(std::string Heat_Prod_Rate_filename){
    
    Array2D net_reaction_rate(n_reaction,points),heat_prod_rate(n_reaction,points);
    vector<double> delta_enthalpy_reaction(n_reaction);
    
        for(size_t j=0;j<points;j++){
            m_flow->setGas( &solution_vector[0],j );
            m_flow->kinetics().getNetRatesOfProgress( &net_reaction_rate(0,j) );
            m_flow->kinetics().getDeltaEnthalpy( &delta_enthalpy_reaction[0] );
            //calculate heat production rate of each reaction
            for(size_t k=0;k<n_reaction;k++)
            {
            heat_prod_rate(k,j)=-net_reaction_rate(k,j)*delta_enthalpy_reaction[k];
            }}
        std::ofstream Heat_prod_rate(Heat_Prod_Rate_filename);
        //print heat production rate to Heat_prod_rate.txt
        for(size_t m=0;m<points;m++){
            //print title
            if(m==0){
                for(size_t n=0;n<n_reaction;n++){
                    //print grid to first column
                    if(n==0){Heat_prod_rate<<std::setw(50)<<std::left<<"Grid[m]"<<"\t";}
                    Heat_prod_rate<<std::setw(50)<<std::left<<reaction_name_vector[n]<<"\t";}
                Heat_prod_rate<<std::endl;}
                
            //print value at each grid point
            for(size_t n=0;n<n_reaction;n++){
                // print grid value to first column
                if(n==0){Heat_prod_rate<<std::setw(50)<<std::left<<std::scientific<<m_flame->domain(1).grid(m);}
                //heat production rate
                Heat_prod_rate<<std::setw(50)<<std::left<<std::scientific<<heat_prod_rate(n,m);
            }
            Heat_prod_rate<<std::endl;
        }
        
    }

void Data_post_analysis::get_temperature_reaction_sensitivity(size_t point_loc){
        //size_t point_loc=50;
        size_t i_Su=m_flame->domain(0).nComponents()+m_flame->domain(1).componentIndex("T")+m_flame->domain(1).nComponents()*point_loc;
        size_t name_id=m_flame->domain(1).componentIndex("T");
        std::cout<<"Get sensitivity of "<<m_flame->domain(1).componentName(name_id)<<" at "<<solution_vector[i_Su]<<std::endl;
        vector<double> dgdx(points*n_variable,0);
        dgdx[i_Su]=1;
        double Su0=solution_vector[i_Su];
        Sens_reaction=solve_adjoint(n_reaction,dgdx);
        for (auto &sens_reaction:Sens_reaction){sens_reaction/=Su0;};
        //for(size_t n=0; n<Sens_reaction.size();n++){if(abs(Sens_reaction[n])>1e-5){std::cout<<std::setw(30)<<std::left<<Sens_reaction[n]<<":"<<std::right<<m_flow->kinetics().reactionString(n)<<"\n";}}
    }

void Data_post_analysis::get_oxygen_sensitivity(size_t point_loc){
        size_t i_Su=m_flame->domain(0).nComponents()+m_flame->domain(1).componentIndex("O2")+m_flame->domain(1).nComponents()*point_loc;
        size_t name_id=m_flame->domain(1).componentIndex("O2");
        std::cout<<"Get sensitivity of "<<m_flame->domain(1).componentName(name_id)<<" at "<<solution_vector[i_Su]<<std::endl;
        vector<double> dgdx(points*n_variable,0);
        dgdx[i_Su]=1;
        double Su0=solution_vector[i_Su];
        Sens_reaction=solve_adjoint(n_reaction,dgdx);
        for (auto &sens_reaction:Sens_reaction){sens_reaction/=Su0;};
      }

void Data_post_analysis::get_sensitivity(size_t point_loc, std::string variable){
        size_t i_Su=m_flame->domain(0).nComponents()+m_flame->domain(1).componentIndex(variable)+m_flame->domain(1).nComponents()*point_loc;
        size_t name_id=m_flame->domain(1).componentIndex(variable);
        vector<double> dgdx(points*n_variable,0);
        std::cout<<"Get sensitivity of "<<m_flame->domain(1).componentName(name_id)<<" at "<<solution_vector[i_Su]<<std::endl;
        dgdx[i_Su]=1;
        double Su0=solution_vector[i_Su];
        Sens_reaction=solve_adjoint(n_reaction,dgdx);
        for (auto &sens_reaction:Sens_reaction){sens_reaction/=Su0;};
      }
    
void Data_post_analysis::get_sensitivity_all(std::string variable){
        Array2D dgdx_all(points*n_variable,points,0.0);
        vector<double> S_var(points);

        for (size_t point_loc=0;point_loc<points;point_loc++){
        size_t i_Su=m_flame->domain(0).nComponents()+m_flame->domain(1).componentIndex(variable)+m_flame->domain(1).nComponents()*point_loc;
        size_t name_id=m_flame->domain(1).componentIndex(variable);
        vector<double> dgdx(points*n_variable,0);
        std::cout<<"Ready to get sensitivity (all) of "<<m_flame->domain(1).componentName(name_id)<<" at "<<solution_vector[i_Su]<<std::endl;
        dgdx_all(i_Su,point_loc)=1;
        S_var[point_loc]=solution_vector[i_Su];}

        solve_adjoint_all(n_reaction,&dgdx_all);

        for (size_t point_loc=0;point_loc<points;point_loc++){
        for (size_t n=0;n<Sens_reaction_all.nRows();n++){Sens_reaction_all.value(n,point_loc)/=S_var[point_loc];};}
    }
    
vector<double> Data_post_analysis::solve_adjoint( size_t n_params, vector<double> dgdx, double(*g_ptr)(), double dp )
    {
        size_t n_vars=points*n_variable;
        vector<double> L(n_vars), fplus(n_vars), fminus(n_vars), dfdp_col(n_vars), dgdp(n_params), L_dfdp(n_params), result(n_params);

        //Array2D dfdp(n_vars, n_params);
        m_flame->solveAdjoint(&dgdx[0],&L[0]);
        std::cout<<"Result from SolveAdjoint is obtained\n";
        //for(auto &e:L){std::cout<<e<<" ";}
        //std::cout<<std::endl;
        clock_t start=clock();
        double g_plus=0, g_minus=0;
        if(dfdp_flag==false){
            dfdp.resize(n_vars,n_params);
        for(size_t n=0;n<n_params;n++) {
            if(n%(n_params/10)==0){std::cout<<"Progress in Perturb: "<<n/(n_params/10)*10<<"%\n";}
            //perturb with dp
            perturb(n,dp);
            g_plus =(g_ptr!=nullptr)? g_ptr():0;
            m_flame->getResidual(0,&fplus[0]);
            //perturb with -dp
            perturb(n,-dp);
            g_minus=(g_ptr!=nullptr)? g_ptr():0;
            m_flame->getResidual(0,&fminus[0]);
            
            //perturb back to 0
            perturb(n,0);
            dgdp[n]=(g_plus-g_minus)/(2*dp);
            for(size_t j=0;j<dfdp.nRows();j++){dfdp(j,n)=(fplus[j]-fminus[j])/(2*dp);}
        }
        /*std::ofstream dfdp_file(folder_name+"dfdp.txt");
        for(size_t l=0;l<n_params;l++){
        for(size_t k=0;k<points;k++){
        for(size_t m=0;m<n_variable;m++){dfdp_file<<std::setw(20)<<std::left<<std::scientific<<dfdp(n_variable*k+m, l);}
        dfdp_file<<std::endl;
        }}*/
            dfdp_flag=true;}
       else {std::cout<<"dfdp has been solved previously\n";}
       for(size_t n=0;n<n_params;n++) {    
            dfdp.getColumn(n,&dfdp_col[0]);
            double sum=0;//obtain inner product between Lambda and dfdp
            for (size_t m=0;m<dfdp_col.size();m++){
                //if(dfdp_col[m]){std::cout<<m<<" "<<dfdp_col[m]<<" "<<L[m]<<std::endl;}
                sum+=dfdp_col[m]*L[m];
            }
            L_dfdp[n]=sum;
            result[n]=dgdp[n]-L_dfdp[n];
        }
        std::cout<<"Time in perturb: "<<double(clock()-start)/CLOCKS_PER_SEC<<"s\n";
        return result;
    }
    
void Data_post_analysis::solve_adjoint_all( size_t n_params, Array2D* dgdx_all, double(*g_ptr)(), double dp )
    {
        size_t n_vars=points*n_variable;
        vector<double> L(n_vars), fplus(n_vars), fminus(n_vars), dfdp_col(n_vars), dgdp(n_params), L_dfdp(n_params), result(n_params);
        Array2D L_all(n_vars,points);
        Sens_reaction_all.resize(n_vars,points,0);
        //Array2D dfdp(n_vars, n_params);
        m_flame->solveAdjoint_all(dgdx_all,&L_all);
        std::cout<<"Result from SolveAdjoint is obtained\n";

        clock_t start=clock();
        double g_plus=0, g_minus=0;
        if(dfdp_flag==false){
            dfdp.resize(n_vars,n_params);
            for(size_t n=0;n<n_params;n++) {
            if(n%(n_params/10)==0){std::cout<<"Progress in Perturb: "<<n/(n_params/10)*10<<"%\n";}
                //perturb with dp
                perturb(n,dp);
                g_plus =(g_ptr!=nullptr)? g_ptr():0;
                m_flame->getResidual(0,&fplus[0]);
                //perturb with -dp
                perturb(n,-dp);
                g_minus=(g_ptr!=nullptr)? g_ptr():0;
                m_flame->getResidual(0,&fminus[0]);

                //perturb back to 0
                perturb(n,0);
                dgdp[n]=(g_plus-g_minus)/(2*dp);
                for(size_t j=0;j<dfdp.nRows();j++){dfdp(j,n)=(fplus[j]-fminus[j])/(2*dp);}
            }
            dfdp_flag=true;}
        else {std::cout<<"dfdp has been solved previously\n";}
        std::cout<<"Time in perturb: "<<double(clock()-start)/CLOCKS_PER_SEC<<"s\n";
        for(size_t j=0;j<points;j++){
            L_all.getColumn(j,&L[0]);

            for(size_t n=0;n<n_params;n++) {
                dfdp.getColumn(n,&dfdp_col[0]);
                double sum=0;//obtain inner product between Lambda and dfdp
                for (size_t m=0;m<dfdp_col.size();m++){
                    sum+=dfdp_col[m]*L[m];
                }
                L_dfdp[n]=sum;
                Sens_reaction_all.value(n,j)=dgdp[n]-L_dfdp[n];
        }}

    }

}
