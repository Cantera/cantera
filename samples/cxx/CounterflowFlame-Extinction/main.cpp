#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/DomainFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/base/Solution.h"
#include "cantera/base/stringUtils.h"
#include "cantera/transport.h"
#include "cantera/onedim.h"
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include "cantera/oneD/Phase_liquid.h"
#include "cantera/oneD/Data_post_analysis.h"
#include "input.dic"
#include "output.h"

using namespace Cantera;
using fmt::print;

std::string construct_stream(const std::map<std::string, double>& fuel_fractions) {
    std::ostringstream stream_builder;
    for (const auto& [species, fraction] : fuel_fractions) 
       { stream_builder << species << ":" << fraction << ",";}
    // Convert to string and remove trailing comma
    std::string stream = stream_builder.str();
    if (!stream.empty()) { stream.pop_back();}
    return stream;}

int extinction(Parameters PA, int solution_status,double increase_SR, std::string *name, std::string *id, int loglevel, std::string file_id_restore="")
{
    clock_t start, start_present, end;
    start=clock();
    try {
         
        std::vector<std::string> fuel_species, ox_species;
        for (const auto& [key, _] : PA.fuel_fractions) {fuel_species.push_back(key);}
        for (const auto& [key, _] : PA.ox_fractions) {ox_species.push_back(key);}
        
        std::string f1 = fuel_species.size() > 0 ? fuel_species[0] : "";
        std::string f2 = fuel_species.size() > 1 ? fuel_species[1] : "";
        std::string f3 = fuel_species.size() > 2 ? fuel_species[2] : "";
        std::string f_stream = construct_stream(PA.fuel_fractions);
        
        std::string ox1 = ox_species.size() > 0 ? ox_species[0] : "";
        std::string ox2 = ox_species.size() > 1 ? ox_species[1] : "";
        std::string ox_stream = construct_stream(PA.ox_fractions);
        std::cout << "f_stream: " << f_stream << " ox_stream: " << ox_stream << std::endl;
        
        auto sol =  newSolution(PA.mech, "gas", "none");
        auto gas = sol->thermo();//solution.h公共函数，返回solution.h中(protected) m_thermo (thermo_phase类型的共享指针)
        double temp_f = PA.T_f, temp_ox = PA.T_ox;//298.15; // Room Temperature
        double SR_ox ; 
        if(solution_status==1) {SR_ox = PA.InitialStrainRate;}// cm/s initial solve without status lablel
        if(solution_status==1 && !file_id_restore.empty()){SR_ox = std::stod(file_id_restore);} //retore from backup id
        if(solution_status==-1){SR_ox = std::stof(*id); file_id_restore="";}// resolve after failure of solving
        
        double pressure = PA.pressure; //atm
        double L = PA.distance;//[m]separation distance 14.5mm
        double u_ox = -SR_ox * L / 4;// Vox = strain rate *L /4 [m/s]

        size_t nsp = gas->nSpecies();//size_t similar to unsigned_int,输出species数量
        vector<double> x_f(nsp, 0.0), y_f(nsp);//本质上依旧是vector, 但用type_def进行替换， 此处生成nsp个元素，均为0. //定义入口质量分数 yin    
        gas->setMixtureFraction(1, f_stream, ox_stream);// =mixFrac*Fuelcomp+(1-mixFrac)*OxComp,
        gas->setState_TP(temp_f, pressure);
        gas->getMoleFractions(x_f.data());//把MoleFraction输出x.data(), 其值为首元素指针
        gas->getMassFractions(&y_f[0]);//把MassFraction输出到&yin[0],依然为yin首元素的地址
        
        double rho_f = gas->density();
        for (int n = 0; n < x_f.size(); n++) {std::cout << x_f[n] << "\t";} std::cout << std::endl;

        vector<double> y_ox(nsp), x_ox(nsp);//定义出口质量分数 //定义出口摩尔分数 新加的
        gas->setMixtureFraction(0, f_stream, ox_stream);// oxidizer gas
        gas->setState_TP(temp_ox, pressure);// 新加的
        gas->getMassFractions(&y_ox[0]);//把平衡时的质量分数设为出口
        gas->getMoleFractions(&x_ox[0]);//新加的
        
        double rho_ox = gas->density();
        std::cout << "rho_f" << rho_f <<  " rho_ox" << rho_ox << std::endl;
        for (int n = 0; n < y_ox.size(); n++) {std::cout << y_ox[n] << "\t";} std::cout << std::endl;
        
        //=============  build each domain ========================

        //-------- step 1: create the flow -------------
        //domain 1: 主要反应区， 需要设置四样：setupGrid, setupTransport, setupKinectics, setupPressure
        gas->setMixtureFraction(0.1,  f_stream, ox_stream);
        gas->equilibrate("HP");//把phase设置到平衡状态，HP为平衡设置关键词
        double T_ad = gas->temperature();
        std::cout<< "Tad = "<< T_ad<<std::endl;
        StFlow flow(gas);
        flow.setAxisymmetricFlow();//m_type=cAxisymmetricStagnationFlow,m_dovisc = true
        // create an initial grid
        int nz = 11;//Grid 数目
        double lz = L;//0.1; 反应区宽度 单位是m
        vector<double> z(nz);//vector 包含grid
        double dz = lz / ((double) (nz - 1));//各个grid的距离=反应区宽度/grid点数
        for (int iz = 0; iz < nz; iz++) {z[iz] = ((double) iz) * dz;}//把grid点输入vector中
        
        flow.setupGrid(nz, &z[0]);//设置初始网格  nz是网格数， &z[0]是网格位置 把&z中的值输入m_z以及m_dz

        // specify the objects to use to compute kinetic rates and
        // transport properties

	std::unique_ptr<Transport> trmix(newTransportMgr("Mix",sol->thermo().get()));//Create a new Transport instance.返回transport object for phase (Transport类指针)
        std::unique_ptr<Transport> trmulti(newTransportMgr("Multi", sol->thermo().get()));//用newTransportMgr返回的结果（Transport类指针）来给指针trimulti赋值

        flow.setTransport(*trmix);
        flow.setKinetics(*sol->kinetics());
        flow.setPressure(pressure);

        //------- step 2: create the inlet_f  -----------------------
        //domain 0 左端入口
        Inlet1D inlet_f;  
        double u_f = -sqrt(rho_ox / rho_f) * u_ox; // u_f is calculated from momentum balance
        double mdot_f = u_f * rho_f;
        std::cout << "m_dot is " << mdot_f << std::endl;
        inlet_f.setMoleFractions(x_f.data());
        inlet_f.setMdot(mdot_f);
        //inlet_f.showMdot();
        inlet_f.setTemperature(temp_f);
        //------- step 3: create the inlet_ox  ---------------------
        //domain 2右端入口
        Inlet1D inlet_ox;
        double mdot_ox = std::abs(u_ox * rho_ox);
        inlet_ox.setMoleFractions(x_ox.data());
        inlet_ox.setMdot(mdot_ox);
        inlet_ox.setTemperature(temp_ox);

        //=================== create the container and insert the domains =====
        //各domain组合并设置值， 需要设置速度，温度，浓度
        std::vector<Domain1D *> domains{&inlet_f, &flow, &inlet_ox};// 从左至右依次排序
	Sim1D flame(domains);
	
        //----------- Supply initial guess----------------------
        //各grid点的位置，按比例计算
        vector<double> locs{0.0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7, 0.8, 0.9, 1.0}, value;
        //对速度进行初步估算
        value = {u_f, u_f / 1.5, u_f / 2, u_f / 2.5, u_f / 4, 0, u_ox / 4, u_ox / 3.5, u_ox / 2, u_ox / 1.5, u_ox};
        flame.setInitialGuess("velocity", locs, value);
        //对温度进行初步估算
        value={temp_f,(0.8*temp_f+0.2*T_ad),(0.5*temp_f+0.5*T_ad), (0.1*temp_f+0.9*T_ad), T_ad, 1.2*T_ad,(0.9*T_ad + 0.1*temp_ox),(0.7*T_ad+0.3*temp_ox), (0.5*T_ad+0.5*temp_ox), (0.2*T_ad+0.8*temp_ox),temp_ox};
        flame.setInitialGuess("T", locs, value);
        for (const auto &T: value) {std::cout << T <<" ";} std::cout<<std::endl; 
        
        //对各物质浓度进行估算 
        for (size_t i = 0; i < nsp; i++) {
            value = {y_f[i], 0.75 * y_f[i], 0.5 * y_f[i], 0.375 * y_f[i], (0.25 * y_f[i]), 0.5 * (y_f[i] + y_ox[i]),(0.25 * y_ox[i]), (0.375 * y_ox[i]), (0.5 * y_ox[i]), (0.75 * y_ox[i]), y_ox[i]};
            flame.setInitialGuess(gas->speciesName(i), locs, value);}
                      
        //设置初始条件
        double moment_bal = rho_f * u_f * u_f / rho_ox / u_ox / u_ox;
        std::cout << "stagnation plane is " << moment_bal << std::endl;
        //inlet_f
        inlet_f.setMoleFractions(x_f.data());//重设入口摩尔分数
        //inlet_ox
        inlet_ox.setMoleFractions(x_ox.data());//重设出口摩尔分数

        int flowdomain = 1;//-1是对所有domain进行refine, 0,1,2 是对domain0,1,2 进行refine
        double ratio = 10,  slope = 0.1, curve = 0.2;//3;//0.1;//0.2;
        flame.setRefineCriteria(flowdomain, ratio, slope, curve);
        flame.setFixedTemperature(1.2*(temp_ox+T_ad));
        
        bool transport=false, soret=false;
        if (transport==true && soret==true) {
                flow.setTransport(*trmulti);
                flow.enableSoret(true);}
        flame.setMaxTimeStepCount(50000);//500
        flame.setMinTimeStep(1e-30);//1e-16 default value 

        //----start loop for searching auto-ignition----up to 1200K
        for (;SR_ox <= 1200;) {
            //flame.showSolution();//show solution of each domain           
            if(solution_status==1 && !file_id_restore.empty()){
                flame.restore("Solutions/solution_"+file_id_restore+"SR.yaml", file_id_restore);
                file_id_restore.clear();
                std::cout<<"--------Restore Solution at "<<std::fixed<<std::setprecision(2)<<SR_ox<<" s^-1 and initial increase of strain rate is "<<increase_SR<<"s^-1--------"<<std::endl;}
            if (solution_status==-1)
            {solution_status=1;
            flame.restore(*name, *id);
            std::cout<<"--------Restore Solution at "<<*id<<" s^-1 and increase of strain rate is "<<increase_SR<<"s^-1--------"<<std::endl;}
             
            //flame.restoreSteadySolution();
            flow.solveEnergyEqn();
            std::cout<<"Try to solve the flame at "<<std::fixed<<std::setprecision(2)<<SR_ox<<" s^-1"<<std::endl;
            
            start_present=clock();
            flame.solve(loglevel, PA.refine_grid);
            end=clock();
            std::cout<<"Time consumed at "<<std::fixed<<std::setprecision(2)<<SR_ox<<" is "<<double(end-start_present)/CLOCKS_PER_SEC/60<<"\n"<<"Total time is "<<double(end-start)/CLOCKS_PER_SEC/60<<std::endl;
            
            OutputFunctions Output(flame, flow, flowdomain, SR_ox, f1, f2, ox1, ox2);          
            double temp_max= Output.Temp_max();
            double diff_T = temp_max-std::max(temp_ox,temp_f);
            Output.printFormattedOutput();  //print to monitor
            Output.yamlsave(name, id,  (diff_T<0.1)); //save solution to yaml file & update name and id [Not update if (diff_T<0.1) is true]
             
            if(PA.save_sol){
            Output.outputcsv();       //save solution to csv file     
            Output.dynamics();}         //print to dynamic file
            
            if (PA.save_eqn==true){
            Output.outputEnergyEqTerms();  //output Energy equation terms
            Output.outputSpeciesEquationTerms();} //output Species equation terms
            
            //check if sensitivity is required
            if(PA.sens_output){Data_post_analysis DA;
            DA.Initialize(flame,flow,sol,std::to_string(temp_ox));}
            
            //check if more calculation is required
            if(increase_SR <= 0.1 && diff_T < 0.1 )
            {std::cout<<"Final Solution is obtained at "<<std::fixed<<std::setprecision(2)<<SR_ox<<std::endl; break;} 
            //adjust the increase of strain rate for next calculation
            if(increase_SR > 0.1 && diff_T < 0.1 ){std::cout<<"Initial Solution is obtained at "<<std::fixed<<std::setprecision(2)<<SR_ox<<std::endl; return 1;}
            //if (diff_T < 1100 && !file_id_restore.empty()) {increase_SR = 5;}
            //update BCs
            SR_ox = SR_ox + increase_SR;
            u_ox = - SR_ox * L/4;
            mdot_ox=std::abs(rho_ox*u_ox);
            u_f = -sqrt(rho_ox / rho_f) * u_ox;
            mdot_f = u_f * rho_f;
            moment_bal = rho_f * u_f * u_f / rho_ox / u_ox / u_ox;
            
            inlet_f.setMdot(mdot_f);
            inlet_f.setTemperature(temp_f);
            inlet_ox.setMdot(mdot_ox);
            inlet_ox.setTemperature(temp_ox);
            
            Output.log(u_ox, increase_SR); //output log file
            } 
       } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << std::endl;
        end=clock();
        std::cout<<"Total time consumed is "<<double(end-start)/CLOCKS_PER_SEC/60<<std::endl;
        return -1;}
        
    end=clock();
    std::cout<<"Total time consumed is "<<double(end-start)/CLOCKS_PER_SEC/60<<std::endl;
    return 0;}

void Get_restore_id_vector(std::string folder_path, std::vector<std::string> &restore_id_vector){
    for(const auto& entry:std::filesystem::directory_iterator(folder_path))
    {if(entry.is_regular_file()&&entry.path().extension()==".yaml")
        {
            std::string filename=entry.path().filename().string();
            std::string id=filename.substr(filename.find("_")+1,filename.find("SR")-filename.find("_")-2);
            restore_id_vector.push_back(id);}
    }
    //order element from low to high
    std::sort(restore_id_vector.begin(),restore_id_vector.end());
    //remove repeated element
    auto it=std::unique(restore_id_vector.begin(),restore_id_vector.end());
    restore_id_vector.erase(it, restore_id_vector.end());
}

void Restore_data_analysis(Parameters PA, std::string solution_path,std::string file_id_restore=""){
    //file_id_restore=-1, restore all solution existed under the folder
    Inlet1D inlet, outlet;
    auto sol = newSolution(PA.mech, "gas", "none");
    auto gas(sol->thermo());
    StFlow flow(gas);
    flow.setAxisymmetricFlow();
    std::unique_ptr<Transport> trmix(newTransportMgr("Mix", sol->thermo().get()));//Create a new Transport instance.返回transport object for phase (Transport类指针)
    flow.setTransport(*trmix);
    flow.setKinetics(*sol->kinetics());
    std::vector<Domain1D *> domains{&inlet, &flow, &outlet};// 从左至右依次排序
    Sim1D flame (domains);

    std::vector<std::string> restore_id_vector;
    std::stringstream id_ss(file_id_restore);
    std::string id;
    if(file_id_restore==""){std::cerr << "No case need to restore" << std::endl;}
    else if(file_id_restore=="-1"){Get_restore_id_vector("Solutions",restore_id_vector);}
    else if(file_id_restore.find(",")!=std::string::npos){while(std::getline(id_ss, id, ',')){id.erase(remove(id.begin(), id.end(),' '),id.end());restore_id_vector.push_back(id);}}
    else {restore_id_vector.push_back(file_id_restore);}
    
    for(size_t k=0;k<restore_id_vector.size();k++){
        std::cout<<k<<" "<<restore_id_vector[k]<<std::endl;
    std::string name="Solutions/solution_"+restore_id_vector[k]+"SR.yaml";
    flame.restore(name,restore_id_vector[k]);
    flame.solve();
    Data_post_analysis DA;
    DA.Initialize(flame,flow, sol, restore_id_vector[k]);
    }}

int main(int argc, char** argv)
{
    Parameters PA;
    int loglevel = 1, status = 1, n = 0;
    std::string sol_name, id, backup_id = PA.backup_id;
    double delta_SR = PA.delta_SR;
    if (PA.mode == "Calculation"){
           while (status!=0){
               status = extinction(PA, status, delta_SR, &sol_name, &id,  loglevel, backup_id);
               backup_id = id;
               delta_SR = delta_SR / std::pow(2,++n); }}
    else if (PA.mode == "Analysis") {Restore_data_analysis(PA,"Solutions", backup_id); }
    
    else {std::cout<<"Please assign a explict mode\n";}
return status;
}
