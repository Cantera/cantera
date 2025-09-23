#include "cantera/base/Parser.h"

#include <fstream>
#include <list>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <cstring>
#include <algorithm>

void Param::print() {
  cout << m_name;
  for (int i = 0; i < m_args.size(); ++i)
     cout << " " << m_args[i];
  cout << endl; 
} 

double Param::readDouble(const int i) {
  double value;
  std::string value_str;
  if ( i >= m_args.size() ) {  
    cout << "---------------------------------------------------------" << endl;
    cout << " Parser error: no argument " << i << " of param " << m_name << endl;
    cout << "---------------------------------------------------------" << endl;
    exit(1);
  } else {
    //std::stringstream stream(m_args[i]);
    value_str = Param::readString(i);
    std::replace( value_str.begin(), value_str.end(), 'd', 'e');
    std::stringstream stream(value_str);
    if ( stream >> value ) {
      return(value); 
    } else {
      cout << "---------------------------------------------------------" << endl;
      cout << " Parser error: argument " << i << " of param " << m_name << " in not a double" << endl;
      cout << "---------------------------------------------------------" << endl;
      exit(1);
    } 
  } 
}

int Param::readInteger(const int i) {
  int value; 
  if ( i >= m_args.size() ) {  
    cout << "---------------------------------------------------------" << endl;
    cout << " Parser error: no argument " << i << " of param " << m_name << endl;
    cout << "---------------------------------------------------------" << endl;
    exit(1);
  } else {
    std::stringstream stream(m_args[i]);
    if ( stream >> value ) {
      return(value); 
    } else {
      cout << "---------------------------------------------------------" << endl;
      cout << " Parser error: argument " << i << " of param " << m_name << " in not an integer" << endl;
      cout << "---------------------------------------------------------" << endl;
      exit(1);
    } 
  } 
}

string Param::readString(const int i) {
  if ( i >= m_args.size() ) {  
    cout << "---------------------------------------------------------" << endl;
    cout << " Parser error: no argument " << i << " of param " << m_name << endl;
    cout << "---------------------------------------------------------" << endl;
    exit(1);
  } else {
    return(m_args[i]);
  } 
}

bool Param::readBool(const int i) {
  if ( i >= m_args.size() ) {  
    cout << "---------------------------------------------------------" << endl;
    cout << " Parser error: no argument " << i << " of param " << m_name << endl;
    cout << "---------------------------------------------------------" << endl;
    exit(1);
  } else {
    if ((m_args[i]=="true")||(m_args[i]=="TRUE")) 
      return(true); 
    else if ((m_args[i]=="false")||(m_args[i]=="FALSE")) 
      return(false); 
    else {
      cout << "---------------------------------------------------------" << endl;
      cout << " Parser error: argument " << i << " of param " << m_name << " in not a boolean" << endl;
      cout << "---------------------------------------------------------" << endl;
      exit(1);
    } 
  } 
}

Param * Parser::getParam(const string& name) {
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if ( p->Name() == name ) {
      p->incrCount();
      return(&(*p));
    } 
  }
  cout << " Param " << name << " is required ! " << endl;
  exit(1);
} 

Param * Parser::getParam(const string& name, const size_t n) {
  size_t count = 0; 
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if (p->Name() == name) {
      if (count == n ) {
        p->incrCount();
        return(&(*p));
      } else {
       ++count; 
      }
    } 
  }
  cout << n <<"th param " << name << " does not exist ! " << endl;
  exit(1);
} 

Param * Parser::getParam(const string& name, const size_t n, const size_t m) {
  size_t global_count = 0;
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if (p->Name() == name && global_count >= n && global_count <= m) {
      p->incrCount(); 
      return(&(*p));
    } 
    global_count = global_count + 1;
  }
  cout << " param " << name << " does not exist for selected mixture ! " << endl;
  exit(1);
} 

// Non-existing method in cantera version
Param * Parser::getParam(const string& name, const size_t n, const size_t m, const size_t l) {
  size_t global_count = 0;
  size_t count = 0; 
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if (p->Name() == name && global_count >=n && global_count <=m) {
      if (count == l ) {
        p->incrCount();
        return(&(*p));
      } else {
       ++count; 
      }
    }
    global_count = global_count + 1; 
  }
  cout << n <<"th param " << name << " does not exist ! " << endl;
  exit(1);
} 


size_t Parser::getParamNumber(const string& name, const size_t n) {
  size_t global_count = 0;
  size_t param_count = 0;
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if (p->Name() == name) {
        if (param_count == n) {
           return(global_count);
	}
        param_count = param_count + 1;
    } 
    global_count = global_count + 1;
  }
  return(string::npos);  
} 

bool Parser::checkParam(const string& name) {
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if ( p->Name() == name ) {
      p->incrCount(); 
      return(true);     
    } 
  }
  return(false); 
} 

// Non-existing method in cantera version
bool Parser::checkParam(const string& name, const size_t n, const size_t m) {
  size_t global_count = 0;
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if ( p->Name() == name && global_count >= n && global_count <= m) {
      p->incrCount(); 
      return(true);     
    }
    global_count = global_count + 1; 
  }
  return(false); 
} 


size_t Parser::nbParamOccurence(const string& name) {
  size_t count = 0;
  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p) {
    if ( p->Name() == name ) {
      ++count;
    }
  }
  return(count);
}

void Parser::parseFile(const string& file) {

// Disclaimer
  //cout << "--------------------------------------------" << endl;
  //cout << " Parsing input file "+ file<< endl;
  //cout << "--------------------------------------------" << endl;

// File exist ?
  ifstream inputfile(file); 
  if (!inputfile.good()) { 
    cout << "------------------------------------------------------" << endl;
    cout << " Parser error: input file "+file+" does not exist !" << endl;
    cout << "------------------------------------------------------" << endl;
    exit(1);
  } 

// Loop on lines
  bool goodParam;
  while (inputfile.good()) {   
    string line;   
    Param buffer_param; 
    goodParam = false;
    getline(inputfile,line);
    if (addParamFromLine(buffer_param,line)) { 
      m_params.push_back(buffer_param);
      m_count += 1; 
      //i_count.push_back(m_count);  
    }
  }

// Print params
//  for (std::list<Param>::iterator p = m_params.begin(); p != m_params.end(); ++p)
//    p->print();

}  

bool Parser::addParamFromLine(Param& param,const string& line) { 
  // look for a hash (Comment operator) and trim  
  size_t hashPos = line.find_first_of("!",0);  
  string trimmed_line = line.substr(0,hashPos);

  // look for first character /= of "="
  string delimiters = " \t=";
  size_t beg = trimmed_line.find_first_not_of(delimiters,0);
   
  // Empty line
  if (beg == string::npos) return(false);

  // look for "=" and store Param name
  size_t delpos = trimmed_line.find_first_of(delimiters,beg);
  param.setName(trimmed_line.substr(beg,delpos-beg));

  // look for parameter arguments, separated by space
  trimmed_line = trimmed_line.substr(delpos+2,trimmed_line.size());

  // Remove initial spaces 
  beg = trimmed_line.find_first_not_of(" ",0);    
  size_t end = trimmed_line.find_last_not_of(" \t");
  trimmed_line = trimmed_line.substr(beg,end);

  // Loop on arguments
  while (!(beg == string::npos)) {
    beg = trimmed_line.find_first_of(" ",0);
    if (!(beg == string::npos)) {
      param.addArgs(trimmed_line.substr(0,beg));
      trimmed_line = trimmed_line.substr(beg,trimmed_line.size());
      size_t space = trimmed_line.find_first_not_of(" ",0);    
      if (!(space == string::npos)) {
        trimmed_line = trimmed_line.substr(space,trimmed_line.size());
      }
    } else {
      param.addArgs(trimmed_line.substr(0,trimmed_line.size()));
    }
  } 
  return(true);
}
