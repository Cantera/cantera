#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/StFlow.h"
namespace Cantera
{
       class OutputFunctions{
       public:
               OutputFunctions(Sim1D &flame, StFlow &flow, int flowdomain,double SR_ox, std::string f1,std::string f2, std::string ox1, std::string ox2){
               m_flame = &flame;
               m_flow = &flow;
               m_flowdomain = flowdomain;
               m_SR_ox = SR_ox;
               m_f1 = f1;
               m_f2 = f2;
               m_ox1 = ox1;
               m_ox2 = ox2;
               prepareVariables();
               
               }
               
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
        
    //else{std::cout<< path << " is already existd\n";}
}

     void printFormattedOutput() {
    // Header
    std::cout << std::setw(12) << "Tf(K)"
              << std::setw(12) << "Tox(K)"
              << std::setw(12) << "Tmax(K)"
              << std::setw(16) << "Tmax-Tox(K)"
              << std::setw(12) << "U_f(m/s)"
              << std::setw(12) << "U_ox(m/s)"
              << std::setw(18) << "rho_f(kg/m^3)"
              << std::setw(18) << "rho_ox(kg/m^3)"
              << std::setw(14) << "MomentBal"
              << std::setw(12) << ("Yf_" + m_f1)
              << std::setw(12) << ("Yf_" + m_f2)
              << std::setw(12) << ("Yox_" + m_ox2)
              << std::setw(12) << ("Yox_" + m_ox1)
              << "\n";

    // Data
    std::cout << std::fixed << std::setprecision(5)
              << std::setw(12) << m_Tvec[0]
              << std::setw(12) << m_Tvec.back()
              << std::setw(12) << m_temp_max
              << std::setw(16) << m_temp_max - m_Tvec.back()
              << std::setw(12) << m_Uvec[0]
              << std::setw(12) << m_Uvec.back()
              << std::setw(18) << m_rhovec[0]
              << std::setw(18) << m_rhovec.back()
              << std::setw(14) << m_moment_bal
              << std::setw(12) << m_Y_f1vec[0]
              << std::setw(12) << m_Y_f2vec[0]
              << std::setw(12) << m_Y_ox2vec.back()
              << std::setw(12) << m_Y_ox1vec.back()
              << "\n";
}           

void prepareVariables() {
    for (int n = 0; n < m_flow->nPoints(); n++) {
        m_Tvec.push_back(m_flame->value(m_flowdomain, m_flow->componentIndex("T"), n));
        m_Uvec.push_back(m_flame->value(m_flowdomain, m_flow->componentIndex("velocity"), n));
        m_Y_f1vec.push_back(m_flame->value(m_flowdomain, m_flow->componentIndex(m_f1), n));
        m_Y_f2vec.push_back(m_flame->value(m_flowdomain, m_flow->componentIndex(m_f2), n));
        m_Y_ox2vec.push_back(m_flame->value(m_flowdomain, m_flow->componentIndex(m_ox2), n));
        m_Y_ox1vec.push_back(m_flame->value(m_flowdomain, m_flow->componentIndex(m_ox1), n));
        m_X_f1vec.push_back(m_flow->Mean_Mw(n) * m_Y_f1vec[n] / m_flow->Mw(m_flow->componentIndex(m_f1) - c_offset_Y));
        m_X_f2vec.push_back(m_flow->Mean_Mw(n) * m_Y_f2vec[n] / m_flow->Mw(m_flow->componentIndex(m_f2) - c_offset_Y));
        m_X_ox2vec.push_back(m_flow->Mean_Mw(n) * m_Y_ox2vec[n] / m_flow->Mw(m_flow->componentIndex(m_ox2) - c_offset_Y));
        m_X_ox1vec.push_back(m_flow->Mean_Mw(n) * m_Y_ox1vec[n] / m_flow->Mw(m_flow->componentIndex(m_ox1) - c_offset_Y));
        m_Condvec.push_back(m_flow->heat_conductivity()[n]);
        m_zvec.push_back(m_flow->grid(n));
        m_rhovec.push_back(m_flow->density(n));
    }
         m_moment_bal =  m_rhovec[0] * m_Uvec[0] * m_Uvec[0] / m_rhovec.back() / m_Uvec.back() / m_Uvec.back();
     m_temp_max = *std::max_element(m_Tvec.begin(), m_Tvec.end());
}

void dynamics()
{
      std::string dynamic_filename = "extinction_"+m_f1+"_"+m_f2+"_massflux_dynamic.csv";
      std::ofstream outfile_dynamic(dynamic_filename, std::ios::out|std::ios::app);
      vector<std::string> title_dynamic={"T_fuel(K)","Tox(K)","Tmax(K)","deltaT(K)","U_f(m/s)","U_ox(m/s)","rho_f(kg/m^3)","rho_ox(kg/m^3)","Y_f"+m_f1,"Y_f"+m_f2,"Y_ox"+m_ox2,"Y_ox"+m_ox1};
      for (size_t n=0;n<title_dynamic.size();n++){outfile_dynamic<<std::setw(18)<<std::left<<title_dynamic[n];}
      outfile_dynamic<<std::endl;
      
      vector<double> output_dynamic={m_Tvec[0], m_Tvec.back(), m_temp_max, m_temp_max- m_Tvec.back(), m_Uvec[0], m_Uvec.back(), m_rhovec[0], m_rhovec.back(), m_Y_f1vec[0], m_Y_f2vec[0], m_Y_ox2vec.back(), m_Y_ox1vec.back()};
        for (size_t m=0; m<output_dynamic.size(); m++){outfile_dynamic<<std::setw(18)<<output_dynamic[m];}
        outfile_dynamic<<std::endl;
        }

void init_outputdynamics(std::ofstream& outfile_dynamic, const std::string f1,const std::string f2,const std::string ox1, const std::string ox2)
{       
        vector<std::string> title_dynamic={"T_fuel(K)","Tox(K)","Tmax(K)","deltaT(K)","U_f(m/s)","U_ox(m/s)","rho_f(kg/m^3)","rho_ox(kg/m^3)","Y_f"+f1,"Y_f"+f2,"Y_ox"+ox2,"Y_ox"+ox1};
        for (size_t n=0;n<title_dynamic.size();n++){outfile_dynamic<<std::setw(18)<<std::left<<title_dynamic[n];}
        outfile_dynamic<<std::endl;
        //outfile_dynamic.close();
        }
        
        
void outputdynamics(std::ofstream& outfile_dynamic,  
                      const std::vector<double>& Tvec, const std::vector<double>& Uvec, 
                      const std::vector<double>& rhovec, const std::vector<double>& Y_f1vec, 
                      const std::vector<double>& Y_f2vec, const std::vector<double>& Y_ox2vec, 
                      const std::vector<double>& Y_ox1vec, double temp_max)
{
        vector<double> output_dynamic={Tvec[0],Tvec.back(),temp_max, temp_max-Tvec.back(), Uvec[0], Uvec.back(), rhovec[0], rhovec.back(), Y_f1vec[0], Y_f2vec[0],Y_ox2vec.back(),Y_ox1vec.back()};
        for (size_t m=0; m<output_dynamic.size(); m++){outfile_dynamic<<std::setw(18)<<output_dynamic[m];}
        outfile_dynamic<<std::endl;
        //outfile_dynamic.close(); 
            }

void outputcsv() {
 check_path("Solution_csv/");
 std::string csv_filename="Solution_csv/extinction_"+m_f1+"_"+m_f2+"_massflux_SR.csv";
 std::string SR_str=std::to_string(m_SR_ox);
 csv_filename.insert(csv_filename.find("SR"),SR_str,0,12);
 std::ofstream csv_outfile(csv_filename, std::ios::trunc);
 //print title to csv file
    vector<std::string>title_outfile = {"Grid(m)", "Temperature(K)", "Uvec(m/s)","rhovec(kg/m^3)","Condvec(W/k-m)","Yf_"+m_f1,"Yf_"+m_f2,"Yo_"+m_ox2, "Yo_"+m_ox1, "Xf_"+m_f1, "Xf_"+m_f2, "Xo_"+m_ox2, "Xo_"+m_ox1};
    for (size_t n=0;n<title_outfile.size();n++){csv_outfile<<std::setw(18)<<std::left<<title_outfile[n];}
            csv_outfile<<std::endl;

    for (size_t n = 0; n < m_flow->nPoints(); n++) {
                vector<double> output={m_flow->grid(n), m_Tvec[n],  m_Uvec[n], m_rhovec[n], m_Condvec[n], m_Y_f1vec[n], m_Y_f2vec[n], m_Y_ox2vec[n], m_Y_ox1vec[n], m_X_f1vec[n], m_X_f2vec[n], m_X_ox2vec[n], m_X_ox1vec[n]};
                
    for (size_t m=0; m<output.size(); m++){csv_outfile<<std::setw(18)<<std::scientific<<std::left<<output[m];}
                csv_outfile<<std::endl;    
               }csv_outfile.close();
}


void outputEnergyEqTerms(){
    check_path("Equation_items/Energy/");
    std::string Eg_filename="Equation_items/Energy/extinction_"+ m_f1 +"_"+ m_f2 +"_massflux_EnegEq_SR.csv";
    std::string SR_str=std::to_string(m_SR_ox);
    Eg_filename.insert(Eg_filename.find("SR"), SR_str, 0, 12);
    outputEnergyEqTerms(Eg_filename, *m_flow);

}

// Function to output Energy Equation terms
void outputEnergyEqTerms( const std::string Eg_filename,  StFlow& flow) {

    // Open the file (trunc mode to overwrite)
    std::ofstream EnegEq_file(Eg_filename, std::ios::trunc);
    if (!EnegEq_file.is_open()) {
        std::cerr << "Error: Unable to open file " << Eg_filename << "\n";
        return;
    }

    // Write the header
    std::vector<std::string> title_outfile = {"Grid(m)", "Term1", "Term2", "Term3", "Term4"};
    for (const auto& title : title_outfile) {
        EnegEq_file << std::setw(18) << std::left << title;
    }
    EnegEq_file << std::endl;

    // Write the data
    for (size_t n = 0; n < flow.nPoints(); n++) {
        std::vector<double> output = {
            flow.grid(n),
            flow.Eneg_terms(0, n),
            flow.Eneg_terms(1, n),
            flow.Eneg_terms(2, n),
            flow.Eneg_terms(3, n)
        };

        for (const auto& value : output) {
            EnegEq_file << std::setw(18) << std::scientific << std::left << value;
        }
        EnegEq_file << std::endl;
    }

    // Close the file
    EnegEq_file.close();
    //std::cout << "Energy equation terms successfully written to " << Eg_filename << "\n";
}


void outputSpeciesEquationTerms (){
    check_path("Equation_items/Species/");
    std::string SR_str=std::to_string(m_SR_ox);
    std::string sp_filename="Equation_items/Species/SpeciesEq.csv";
    outputSpeciesEquationTerms(sp_filename,  SR_str, *m_flow);
}

void outputSpeciesEquationTerms(const std::string& sp_filename, const std::string& SR_str, 
                                 StFlow& flow) {
    // Extract base name and extension
    size_t dot_pos = sp_filename.find_last_of(".");
    std::string base_name = sp_filename.substr(0, dot_pos);
    std::string ext_name = sp_filename.substr(dot_pos);

    // Loop over all species (components)
    for (size_t n = 0; n < flow.nComponents() - c_offset_Y; n++) {
        // Construct file name for each species
        std::string species_filename = base_name + "_" + SR_str + "_" + flow.componentName(n + c_offset_Y) + ext_name;

        // Open the file
        std::ofstream SpEq_file(species_filename);
        if (!SpEq_file.is_open()) {
            std::cerr << "Error: Unable to open file " << species_filename << "\n";
            continue;
        }

        // Write header
        std::vector<std::string> title_outfile = {"Grid(m)", "Convection", "Diffusion", "Generation"};
        for (const auto& title : title_outfile) {
            SpEq_file << std::setw(18) << std::left << title;
        }
        SpEq_file << std::endl;

        // Write data for each grid point
        for (size_t j = 0; j < flow.nPoints(); j++) {
            std::vector<double> output = {
                flow.grid(j),
                flow.Species_terms(0, j, n),
                flow.Species_terms(1, j, n),
                flow.Species_terms(2, j, n)
            };

            for (const auto& value : output) {
                SpEq_file << std::setw(18) << std::scientific << std::left << value;
            }
            SpEq_file << std::endl;
        }

        // Close the file
        SpEq_file.close();}
    }
   
    
void yamlsave(std::string* name, std::string* id, bool updateid = false) {
    // Ensure the base path exists
    check_path("Solutions/");

    // Construct the solution file name
    std::string solutionname = "Solutions/solution_SR.yaml"; // Example: "Solutions/solution_SR.yaml"
    std::string solution_SR = std::to_string(m_SR_ox);
    solutionname.insert(solutionname.find("SR"), solution_SR, 0, 12);
    // Save the solution
    m_flame->save(solutionname, solution_SR, "initial solution", true); // overwrite = true
    // Update the output parameters
    if (!updateid){
    *name = solutionname;
    *id = solution_SR;}
}

//void yamlsave(std::string* name, std::string* id) { yamlsave(name, id, false);}

void log(double u_ox, double increase_SR){
            std::string log="log.txt";
            std::ofstream outfile_log(log, std::ios::out|std::ios::app);
            outfile_log << u_ox<<"\t"<<increase_SR<<"\n";
            outfile_log.close();
}
   double Temp_max(){return m_temp_max;}
   
   
   protected:
        Sim1D *m_flame;
        StFlow *m_flow;
        int m_flowdomain;
        double m_moment_bal, m_temp_max, m_SR_ox;
        vector<double> m_zvec, m_Tvec, m_Uvec, m_rhovec, m_Condvec, m_Y_f1vec, m_Y_f2vec, m_Y_ox2vec, m_Y_ox1vec, m_X_f1vec, m_X_f2vec, m_X_ox2vec, m_X_ox1vec;
        std::string m_f1, m_f2, m_ox1, m_ox2;
   };
}
