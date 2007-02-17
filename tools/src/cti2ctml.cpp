/**
 *  @file example2.cpp
 *
 */

//  Example 
//
//  Read a mechanism and a thermodynamics file for the 
//  class IdealSolidSolnPhase in order to test that it's
//  working correctly
//

#include <iostream>
#include <string>
#include <vector>

#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"

using namespace Cantera;
using namespace std;

#ifdef DEBUG_HKM
int iDebug_HKM = 0;
#endif

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
static void printUsage()
{
    cout << "ctitoxml [-h] infile.cti" << endl;
    cout << "    Translates a cti formated file to an xml file" << endl;
    cout << "    The xml file will be named infile.xml" << endl;
}






int main(int argc, char** argv) {
    string infile;  
    // look for command-line options
    if (argc > 1) {
      string tok;
      for (int j = 1; j < argc; j++) {
	tok = string(argv[j]);
	if (tok[0] == '-') {
	  int nopt = static_cast<int>(tok.size());
	  for (int n = 1; n < nopt; n++) {
	    if (tok[n] == 'h') {
	      printUsage();
	      exit(0);
	    } else {
	      printUsage();
	      exit(1);
	    }
	  }
	} else if (infile == "") {
	  infile = tok;
	}
	else {
	  printUsage();
	  exit(1);
	}
      }
    }
    if (infile == "") {
	  printUsage();
	  exit(1);      
    }
 
    try {
      XML_Node *xc = new XML_Node();
      string path = findInputFile(infile);
      ctml::get_CTML_Tree(xc, path, 0); 
      XML_Node *xd = new XML_Node();
      xc->copy(xd);
      ofstream tout;
      tout.open("testdest.xml");
      xc->write(tout);
      tout.close();
      tout.open("testdest2.xml");
      xd->write(tout);
      tout.close();

    }
    catch (CanteraError) {
      showErrors(cout);
    }

    return 0;
}
/***********************************************************/
