#include "cantera/base/ctexceptions.h"

void demoprog()
{
    // Calls Cantera
}

int main(int argc, char** argv)
{
    try {
        demoprog();
    } catch (Cantera::CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
    return 0;
}
