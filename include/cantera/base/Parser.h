// Small parser for driving CANTERA C++ program
#include <iostream>
#include <list>
#include <vector>
using std::cout;
using std::endl;
using std::string;
#include <fstream>
using std::ifstream;
#include <cstring>

class Param {
private:
    size_t m_count;
    string m_name;
    std::vector<string> m_args;
public:
    Param() {
      m_count = 0;
    }
    Param(const string& name) {
      m_count = 0;
      setName(name);
    }
    void setName(const string& name) {
      m_name = name;
    }
    string Name() {
      return m_name;
    }
    void print(); 
    void addArgs(const string& arg) {
      m_args.push_back(arg);
    } 
    size_t nbArgs() {
      return(m_args.size());
    }
    void incrCount() {
      ++m_count;
    } 
    double readDouble(const int i);    
    int readInteger(const int i);    
    bool readBool(const int i);    
    string readString(const int i);    
};

class Parser {
private:
  std::list<Param> m_params; 
  //std::list<size_t> i_count;
  size_t m_count;
public:
  Parser() {
    m_count = 0;
  }  
  void parseFile(const string& file);
  bool addParamFromLine(Param& param,const string& line);
  Param * getParam(const string& name);
  Param * getParam(const string& name, const size_t n);
  Param * getParam(const string& name, const size_t n, const size_t m);
  Param * getParam(const string& name, const size_t n, const size_t m, const size_t l);
  size_t getParamNumber(const string& name, const size_t n);
  bool checkParam(const string& name);
  bool checkParam(const string& name, const size_t n, const size_t m);
  size_t nbParamOccurence(const string& name);
};
