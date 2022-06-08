/*!
 * @file demprog.cpp
 *
 * Test file demonstrating shared library support - kept in a separate folder to
 * enable testing with different Makefile's
 *
 * Keywords: combustion
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/Solution.h"
#include <iostream>

using namespace Cantera;

void demoprog()
{
    auto sol = newSolution("h2o2.yaml", "ohmech");
}

int main(int argc, char** argv)
{
    try {
        std::cout << "running demoprog..." << std::endl;
        std::cout.flush();
        demoprog();
        std::cout << "finished!" << std::endl;
    } catch (std::exception& err) {
        std::cout << "oops..." << std::endl;
        std::cout << err.what() << std::endl;
    }
}
