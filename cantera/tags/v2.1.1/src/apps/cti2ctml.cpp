/**
 *  @file cti2ctml.cpp
 */

#include "cantera/base/ct_defs.h"
#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"

using namespace Cantera;
using namespace std;


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
static void printUsage()
{
    cout << "cti2ctml [-h] infile.cti" << endl;
    cout << "    Translates a cti formated file to an xml file" << endl;
    cout << "    The xml file will be named ./basename(infile).xml" << endl;
    cout << "     - It will always be written to the current directory" << endl;
}



/*****************************************************************/


int main(int argc, char** argv)
{
    std::string infile;
    // look for command-line options
    if (argc > 1) {
        std::string tok;
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
            } else {
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
        XML_Node* xc = new XML_Node();
        std::string path = findInputFile(infile);
        ctml::get_CTML_Tree(xc, path, 0);
        //XML_Node *xd = new XML_Node();
        //xc->copy(xd);

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }

    return 0;
}
/***********************************************************/
