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

using namespace Cantera;
using fmt::print;

void check_path(std::string path)
{
    struct stat st;
    if(stat(path.c_str(), &st) == -1){//check if path is already existed
    
       size_t pos=0;   
       while((pos = path.find_first_of('/', pos+1))!=std::string::npos){
             std::string subDir = path.substr(0 ,pos);
             if(mkdir(subDir.c_str(), 0700) && errno != EEXIST){break;}}         
       if(mkdir(path.c_str(), 0700) && errno != EEXIST){std::cout<<path<< " can not created\n";}
       else{std::cerr<< path << " has been created \n";;}}
        
    else{std::cout<< path << " is already existd\n";}
}

int auto_ignition( Parameters PA, int solution_status,double increase_temp, std::string *name, std::string *id, int loglevel, std::string file_id_restore="")
{
    clock_t start, start_present, end;
    start=clock();
    try {
        auto sol = newSolution(PA.mech,"gas", "mixture-averaged");
        auto gas = sol->thermo();//solution.h公共函数，返回solution.h中(protected) m_thermo (thermo_phase类型的共享指针)
        double temp_f = 330; //315; // K
        double temp_ox;
        if(solution_status==1) {temp_ox = PA.InitialTemperature;}//480;
        if(solution_status==1 && !file_id_restore.empty()){temp_ox=std::stod(file_id_restore);}
        if(solution_status==-1){temp_ox=std::stof(*id);file_id_restore="";}
        double pressure = PA.pressure; //atm
        double L = PA.distance;
        double u_ox = PA.u_ox;

        size_t nsp = gas->nSpecies();//size_t similar to unsigned_int,输出species数量
        vector<double> x_f(nsp, 0.0);//本质上依旧是vector, 但用type_def进行替换， 此处生成nsp个元素，均为0.
            
        gas->setMixtureFraction(1, "C2H5OH:1e-22, N2:0.6,NC7H16:0.4","O2:0.21,N2:0.79");//"C2H5OH:0.65,N2:0.3,H2O:0.05",
        gas->setState_TP(temp_f, pressure);
        gas->getMoleFractions(x_f.data());//把MoleFraction输出x.data(), 其值为首元素指针
        for (int n = 0; n < x_f.size(); n++) {std::cout << x_f[n] << "\t";}
        std::cout << std::endl;

        double rho_f = gas->density();
        std::cout << "rho_f" << rho_f << std::endl;

        vector<double> y_f(nsp);//定义入口质量分数 yin
        gas->getMassFractions(&y_f[0]);//把MassFraction输出到&yin[0],依然为yin首元素的地址

        //gas->equilibrate("HP");//把phase设置到平衡状态，HP为平衡设置关键词
        vector<double> y_ox(nsp);//定义出口质量分数
        vector<double> x_ox(nsp);//定义出口摩尔分数 新加的

        gas->setMixtureFraction(0, "C2H5OH:0.7,N2:0.3", "O2:0.21,N2:0.79");
        gas->setState_TP(temp_ox, pressure);// 新加的
        gas->getMassFractions(&y_ox[0]);//把平衡时的质量分数设为出口
        gas->getMoleFractions(&x_ox[0]);//新加的
        double rho_ox = gas->density();
        std::cout << "rho_ox" << rho_ox << std::endl;
            
        for (int n = 0; n < y_ox.size(); n++) {std::cout << y_ox[n] << "\t";}
        std::cout << std::endl;
        
        double u_f = -sqrt(rho_ox / rho_f) * u_ox / 2;//default and valid at 800K
        //=============  build each domain ========================

        //-------- step 1: create the flow -------------
        //domain 1: 主要反应区， 需要设置四样：setupGrid, setupTransport, setupKinectics, setupPressure

        StFlow flow(gas);
        flow.setAxisymmetricFlow();//m_type=cAxisymmetricStagnationFlow,m_dovisc = true
        // create an initial grid
        int nz = 11;//Grid 数目
        double lz = 0.0105;//0.1; 反应区宽度 单位是m
        vector<double> z(nz);//vector 包含grid
        double dz = lz / ((double) (nz - 1));//各个grid的距离=反应区宽度/grid点数
        for (int iz = 0; iz < nz; iz++) {z[iz] = ((double) iz) * dz;}//把grid点输入vector中

        flow.setupGrid(nz, &z[0]);//设置初始网格  nz是网格数， &z[0]是网格位置 把&z中的值输入m_z以及m_dz

        // specify the objects to use to compute kinetic rates and
        // transport properties

        std::unique_ptr<Transport> trmix(newTransportMgr("Mix",sol->thermo().get()));//Create a new Transport instance.返回transport object for phase (Transport类指针)
        std::unique_ptr<Transport> trmulti(
                    newTransportMgr("Multi", sol->thermo().get()));//用newTransportMgr返回的结果（Transport类指针）来给指针trimulti赋值
        flow.setTransport(*trmix);
        flow.setKinetics(*sol->kinetics());
        flow.setPressure(pressure);
        //flow.Enable_GPU(true);//enbale GPU for multicomponents transport

        //------- step 2: create the inlet_f  -----------------------
        //domain 0 左端入口
        Inlet1D_new inlet_f;
        double mdot_f = u_f * rho_f;
        std::cout << "m_dot is " << mdot_f << std::endl;
        inlet_f.setMoleFractions(x_f.data());
        inlet_f.setMdot(mdot_f);
        inlet_f.showMdot();
        inlet_f.setTemperature(temp_f);

        Phase_liquid pool;//add liquid pool to boundary
        pool.setMoleFraction(PA.pool_molar_fraction);
        inlet_f.addLiquidBcs(pool);

        print("Inlet_f is set \n");

        //------- step 3: create the inlet_ox  ---------------------
        //domain 2右端入口
        Inlet1D_new inlet_ox;
        double mdot_ox = std::abs(u_ox * rho_ox);
        inlet_ox.setMoleFractions(x_ox.data());
        inlet_ox.setMdot(mdot_ox);
        inlet_ox.setTemperature(temp_ox);

        //=================== create the container and insert the domains =====
        //各domain组合并设置值， 需要设置速度，温度，浓度
        std::vector<Domain1D *> domains{&inlet_f, &flow, &inlet_ox};// 从左至右依次排序
        Sim1D flame(domains);
        //----------- Supply initial guess----------------------

        vector<double> locs{0.0, 0.01, 0.015,0.02,0.025,0.03,0.035,0.15, 0.2, 0.6, 1.0};//各grid点的位置，按比例计算
        vector<double> value;

        value = {u_f, u_f / 1.5, u_f / 2, u_f / 2.5, u_f / 4, 0, u_ox / 4, u_ox / 3.5, u_ox / 2, u_ox / 1.5, u_ox};
        flame.setInitialGuess("velocity", locs, value);//对速度进行初步估算

        value={temp_f,(0.8*temp_f+0.2*temp_ox),(0.6*temp_f+0.4*temp_ox),(0.4*temp_f+0.6*temp_ox),(0.2*temp_f+0.8*temp_ox),temp_ox,temp_ox,temp_ox,temp_ox,temp_ox,temp_ox};

        for (auto &T: value) {std::cout << T << std::endl;}
            
        flame.setInitialGuess("T", locs, value);//对温度进行初步估算
        for (size_t i = 0; i < nsp; i++) {
            value = {y_f[i], 0.75 * y_f[i], 0.5 * y_f[i], 0.375 * y_f[i], (0.25 * y_f[i]), 0.5 * (y_f[i] + y_ox[i]),
                         (0.25 * y_ox[i]), (0.375 * y_ox[i]), (0.5 * y_ox[i]), (0.75 * y_ox[i]), y_ox[i]};
            flame.setInitialGuess(gas->speciesName(i), locs, value);//对各物质浓度进行估算
        }
        //设置初始条件
        std::cout << "stagnation plane is " << rho_f * u_f * u_f / rho_ox / u_ox / u_ox << std::endl;
        //inlet_f
        inlet_f.setMoleFractions(x_f.data());//重设入口摩尔分数
        //inlet_ox
        inlet_ox.setMoleFractions(x_ox.data());//重设出口摩尔分数

        bool transport=false, soret=false;
        
        std::string dynamic_filename="autoignition_C2H5OH_Hep_massflux_dynamic.csv";//for dynamic file
        std::ofstream outfile_dynamic(dynamic_filename, std::ios::out|std::ios::app);
        vector<std::string> title_dynamic={"T_fuel(K)","Tox(K)","Tmax(K)","deltaT(K)","U_f(m/s)","U_ox(m/s)","rho_f(kg/m^3)","rho_ox(kg/m^3)","Y_f(C2H5OH)","Y_f(C7H16)","Yf_lout(C2H5OH)","Yf_lout(C7H16)","Y_ox(N2)","Y_ox(O2)"};
        for (size_t n=0;n<title_dynamic.size();n++){outfile_dynamic<<std::setw(18)<<std::left<<title_dynamic[n];}
        outfile_dynamic<<std::endl;
        outfile_dynamic.close();

        //----start loop for searching auto-ignition----up to 1200K
        for (;temp_ox<=1300;) {
            inlet_f.setMdot(mdot_f);
            inlet_f.setTemperature(temp_f);
            inlet_ox.setMdot(mdot_ox);
            inlet_ox.setTemperature(temp_ox);
            //flame.showSolution();//show solution of each domain

            int flowdomain = 1;//-1是对所有domain进行refine, 0,1,2 是对domain0,1,2 进行refine
            double ratio = 3;//2;//3; (>2)
            double slope = 0.1;//0.025;//0.1;
            double curve = 0.2;//0.05;//0.2;

            flame.setRefineCriteria(flowdomain, ratio, slope, curve);
            // Solve flame
            flow.solveEnergyEqn();
            if (transport==true && soret==true)
            {
                flow.setTransport(*trmulti);
                flow.enableSoret(true);
            }
            flame.setMaxTimeStepCount(5000);//1000//500
            flame.setMinTimeStep(1e-50);//1e-50//1e-16 default value
            if(solution_status==1 && !file_id_restore.empty())
            {
                flame.restore("Solutions/solution_"+file_id_restore+"K.yaml", file_id_restore);
                file_id_restore.clear();
                std::cout<<"--------Restore Solution at "<<std::fixed<<std::setprecision(6)<<temp_ox<<"K and initial increase of temperature is "<<increase_temp<<"K--------"<<std::endl;}
            if (solution_status==-1)
            {solution_status=1;
            flame.restore(*name, *id);
             std::cout<<"--------Restore Solution at "<<*id<<"K and increase of temperature is "<<increase_temp<<"K--------"<<std::endl;}
            start_present=clock();
            flame.solve(loglevel, PA.refine_grid);
            //flame.writeStats(1);
            //flame.showSolution();
            
            end=clock();
            std::cout<<"Time consumed at "<<std::fixed<<std::setprecision(4)<<temp_ox<<" is "<<double(end-start_present)/CLOCKS_PER_SEC/60<<"\n"<<"Total time is "<<double(end-start)/CLOCKS_PER_SEC/60<<std::endl;
            
            //prepare varibles for quick output
            vector<double> zvec,Tvec,Uvec,rhovec,Condvec,Y_C2H5OHvec,Y_C7H16vec,Y_N2vec,Y_O2vec,X_C2H5OHvec,X_C7H16vec,X_N2vec,X_O2vec;
            for (size_t n = 0; n < flow.nPoints(); n++) {
                Tvec.push_back(flame.value(flowdomain,flow.componentIndex("T"),n));
                Uvec.push_back(flame.value(flowdomain,flow.componentIndex("velocity"),n));
                Y_C2H5OHvec.push_back(flame.value(flowdomain,flow.componentIndex("C2H5OH"),n));
                Y_C7H16vec.push_back(flame.value(flowdomain,flow.componentIndex("NC7H16"),n));
                Y_N2vec.push_back(flame.value(flowdomain,flow.componentIndex("N2"),n));
                Y_O2vec.push_back(flame.value(flowdomain,flow.componentIndex("O2"),n));
                X_C2H5OHvec.push_back(flow.Mean_Mw(n)*Y_C2H5OHvec[n]/flow.Mw(flow.componentIndex("C2H5OH")-c_offset_Y));
                X_C7H16vec.push_back(flow.Mean_Mw(n)*Y_C7H16vec[n]/flow.Mw(flow.componentIndex("NC7H16")-c_offset_Y));
                X_N2vec.push_back(flow.Mean_Mw(n)*Y_N2vec[n]/flow.Mw(flow.componentIndex("N2")-c_offset_Y));
                X_O2vec.push_back(flow.Mean_Mw(n)*Y_O2vec[n]/flow.Mw(flow.componentIndex("O2")-c_offset_Y));
                Condvec.push_back(flow.heat_conductivity()[n]);
                zvec.push_back(flow.grid(n));
                rhovec.push_back(flow.density(n));
            }
            double delta_H=pool.Delta_H_mass(Tvec[0]);
            double temp_max=*std::max_element(Tvec.begin(), Tvec.end());
            double diff_T=temp_max-std::max(temp_ox,temp_f);
            
             //print to monitor
             size_t n=(flow.nPoints()-1);
             print("\n{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\t{:12s}\n",
                  "Tf (K)","Tox (K)", "Tmax (K)","Tox_Tmax(K)", "U_f(m/s)","U_ox(m/s)","rho_f(kg/m^3)","rho_ox(kg/m^3)","Y_f(C2H5OH)","Y_f(C7H16)", "Y_ox(N2)", "Y_ox(O2)");
             print("{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.5f}\t{:12.6f}\t{:12.6f}\t{:12.6f}\t{:12.6f}\t{:12.6f}\t{:12.6f}\t{:12.6f}\t{:12.6f}\n",
                Tvec[0],Tvec[n], temp_max, temp_max-Tvec[n], Uvec[0], Uvec[n], rhovec[0], rhovec[n], Y_C2H5OHvec[0], Y_C7H16vec[0], Y_N2vec[n], Y_O2vec[n]);
             printf("Max temperature is %12.6f and T_max-T_ox is %12.6f\n", temp_max, temp_max-Tvec[n]);
                
            if(PA.save_sol){
            //save solution to yaml file
            check_path("Solutions/");
            std::string solutionname="Solutions/solution_K.yaml";
            std::string solution_temp=std::to_string(temp_ox);
            solutionname.insert(solutionname.find("K"),solution_temp,0,12);
            flame.save(solutionname,solution_temp, "initial solution", true);
            *name=solutionname;
            *id=solution_temp;
            
            //save solution to csv file
            check_path("Solution_csv/");
            std::string filename="Solution_csv/autoignition_C2H5OH_Hep_massflux_K.csv";
            std::string temp_str=std::to_string(temp_ox);
            filename.insert(filename.find("K"),temp_str,0,12);
            std::ofstream outfile(filename, std::ios::trunc);
            
            //print title to csv file
            vector<std::string>title_outfile={"Grid(m)", "Temperature(K)", "Uvec(m/s)","rhovec(kg/m^3)","Condvec(W/k-m)","delta_H(J/kg)","Y(C2H5OH)","Y(C7H16)","Y(N2)", "Y(O2)", "X(C2H5OH)", "X(C7H16)", "X(N2)", "X(O2)", "X_lout(C2H5OH)", "X_lout(C7H16)"};
            for (size_t n=0;n<title_outfile.size();n++){outfile<<std::setw(18)<<std::left<<title_outfile[n];}
            outfile<<std::endl;

            for (size_t n = 0; n < flow.nPoints(); n++) {
                vector<double> output={flow.grid(n), Tvec[n],  Uvec[n], rhovec[n], Condvec[n], delta_H, Y_C2H5OHvec[n], Y_C7H16vec[n], Y_N2vec[n], Y_O2vec[n], X_C2H5OHvec[n], X_C7H16vec[n], X_N2vec[n], X_O2vec[n], pool.MoleFraction_outlet("C2H5OH"), pool.MoleFraction_outlet("NC7H16")};
                for (size_t m=0;m<output.size();m++){outfile<<std::setw(18)<<std::scientific<<std::left<<output[m];}
                outfile<<std::endl;    
               }outfile.close();
               
             //print to dynamic file
             size_t n=(flow.nPoints()-1);
             outfile_dynamic.open("autoignition_C2H5OH_Hep_massflux_dynamic.csv",std::ios::app);
             vector<double> output_dynamic={Tvec[0],Tvec[n],temp_max,temp_max-Tvec[n], Uvec[0], Uvec[n], rhovec[0], rhovec[n], Y_C2H5OHvec[0], Y_C7H16vec[0],pool.MoleFraction_outlet("C2H5OH"), pool.MoleFraction_outlet("NC7H16"),Y_N2vec[n],Y_O2vec[n]};
             for (size_t m=0; m<output_dynamic.size(); m++){outfile_dynamic<<std::setw(18)<<output_dynamic[m];}
             outfile_dynamic<<std::endl;
             outfile_dynamic.close(); 
            }
            
            if (PA.save_eqn==true){
            //output Energy equation terms
            std::string temp_folder="Equation_items/Energy/";
            std::string temp_str=std::to_string(temp_ox);           
            std::string temp_filename=temp_folder+"autoignition_C2H5OH_Hep_massflux_EnegEq_"+temp_str+"K.csv";
            check_path(temp_folder);
            std::ofstream EnegEq_file(temp_filename);

            std::vector<std::string> title_outfile={"Grid(m)", "Term1", "Term2",  "Term3", "Term4"};
            for(size_t n=0;n<title_outfile.size();n++){EnegEq_file<<std::setw(18)<<std::left<<title_outfile[n];}
            EnegEq_file<<std::endl;

            for (size_t n = 0; n < flow.nPoints(); n++) {
                vector<double> output={flow.grid(n), flow.Eneg_terms(0,n),flow.Eneg_terms(1,n),flow.Eneg_terms(2,n),flow.Eneg_terms(3,n)};
                for(size_t m=0;m<output.size();m++){EnegEq_file<<std::setw(18)<<std::scientific<<std::left<<output[m];}
                EnegEq_file<<std::endl;
            }EnegEq_file.close();
            
            //output Species equation terms
            std::string sp_folder="Equation_items/Species/";
            std::string sp_filename=sp_folder+temp_str+"K/SpeciesEq.csv";
            check_path(sp_folder+temp_str+"K");
            
            size_t dot_pos=sp_filename.find_last_of(".");
            std::string base_name=sp_filename.substr(0,dot_pos);
            std::string ext_name=sp_filename.substr(dot_pos);
            
            for (size_t n = 0; n < flow.nComponents()-c_offset_Y; n++){
            std::ofstream SpEq_file(base_name + "_" + temp_str + "K_" + flow.componentName(n+c_offset_Y) + ext_name);
            std::vector<std::string> title_outfile={"Grid(m)", "Convection", "Diffusion", "Generation"};
            for(size_t n=0;n<title_outfile.size();n++){SpEq_file<<std::setw(18)<<std::left<<title_outfile[n];}
            SpEq_file<<std::endl;

            for (size_t j = 0; j < flow.nPoints(); j++) {
                vector<double> output={flow.grid(j), flow.Species_terms(0, j, n), flow.Species_terms(1, j, n), flow.Species_terms(2, j, n)};
                for(size_t m=0;m<output.size();m++){SpEq_file<<std::setw(18)<<std::scientific<<std::left<<output[m];}
                SpEq_file<<std::endl;
            }SpEq_file.close();}
            }
            
            //check if sensitivity is required
            if(PA.sens_output){Data_post_analysis DA;
            DA.Initialize(flame,flow,inlet_f,inlet_ox,sol,std::to_string(temp_ox));}
            
            //check if more calculation is required
            if(diff_T>100 )
            {std::cout<<"Final Solution is obtained at "<<std::fixed<<std::setprecision(4)<<temp_ox<<std::endl;
                break;} 
            //adjust the increase of temperature for next calculation
            if (diff_T>0.0001&&diff_T<1) {increase_temp=1<increase_temp? 1:increase_temp;}
            if (diff_T>=1) {increase_temp=0.1<increase_temp? 0.1:increase_temp;}
            if (diff_T>=10) {increase_temp=0.01<increase_temp? 0.01:increase_temp;}
            
            rho_ox=rho_ox*temp_ox/(temp_ox+increase_temp);
            mdot_ox=std::abs(rho_ox*u_ox);
            temp_ox=temp_ox+increase_temp;
            flame.restoreSteadySolution();
            std::cout<<"Try to solve the flame at "<<std::fixed<<std::setprecision(4)<<temp_ox<<" K"<<std::endl;
            
            std::string log="log.txt";//for dynamic file
            std::ofstream outfile_log(log, std::ios::out|std::ios::app);
            outfile_log << temp_ox<<"\t"<<increase_temp<<"\n";
            outfile_log.close();
        }
} catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << std::endl;
        end=clock();
        std::cout<<"Total time consumed is "<<double(end-start)/CLOCKS_PER_SEC/60<<std::endl;
        return -1;
    }
    end=clock();
    std::cout<<"Total time consumed is "<<double(end-start)/CLOCKS_PER_SEC/60<<std::endl;
    return 0;
}

void Get_restore_id_vector(std::string folder_path, std::vector<std::string> &restore_id_vector){
    for(const auto& entry:std::filesystem::directory_iterator(folder_path))
    {if(entry.is_regular_file()&&entry.path().extension()==".yaml")
        {
            std::string filename=entry.path().filename().string();
            std::string id=filename.substr(filename.find("_")+1,filename.find("K")-filename.find("_")-1);
            restore_id_vector.push_back(id);}
    }
    //order element from low to high
    std::sort(restore_id_vector.begin(),restore_id_vector.end(),[](const std::string& a, const std::string& b){return std::stod(a)<std::stod(b);});
    //remove repeated element
    auto it=std::unique(restore_id_vector.begin(),restore_id_vector.end());
    restore_id_vector.erase(it, restore_id_vector.end());
}

void Restore_data_analysis(Parameters PA, std::string solution_path,std::string file_id_restore=""){
    //file_id_restore=-1, restore all solution existed under the folder
    Inlet1D_new inlet, outlet;
    Phase_liquid pool;//add liquid pool to boundary
    pool.setMoleFraction(PA.pool_molar_fraction);
    inlet.addLiquidBcs(pool);
    auto sol = newSolution(PA.mech, "gas", "none");
    auto gas(sol->thermo());
    StFlow flow(gas);
    flow.setAxisymmetricFlow();
    std::unique_ptr<Transport> trmix(newTransportMgr("Mix",sol->thermo().get()));//Create a new Transport instance.返回transport object for phase (Transport类指针)
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
        std::cout<<k+1<<"/"<<restore_id_vector.size()<<" "<<restore_id_vector[k]<<std::endl;
    std::string name="Solutions/solution_"+restore_id_vector[k]+"K.yaml";
    flame.restore(name,restore_id_vector[k]);
    //flow.Enable_GPU(true);//enbale GPU for multicomponents transport
    //flame.solve();
    Data_post_analysis DA;
    DA.Initialize(flame,flow,inlet, outlet, sol, restore_id_vector[k]);
    }
}

int main(int argc, char** argv)
{   
    Parameters PA;
    
    int loglevel = 1, status=1, n=0;
    std::string sol_name, id, backup_id=PA.backup_id;
    double delta_T=PA.delta_T;
    if (PA.mode == "Calculation"){
        while (status!=0){
            status=auto_ignition( PA, status, delta_T, &sol_name, &id, loglevel, backup_id);
            n++;
            delta_T=delta_T/std::pow(1,n);}}

    else if (PA.mode == "Analysis"){
        Restore_data_analysis(PA, "Solutions", backup_id);}

    else {std::cout<<"Please assign a explict mode\n";}
return status;
}
