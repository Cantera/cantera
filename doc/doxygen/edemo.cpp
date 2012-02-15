#include <cantera/base/ctexceptions.h>

void demoprog()
{
    // Calls Cantera
}

int main(int argc, char** argv)
{
    try {
        demoprog();
    } catch (Cantera::CanteraError) {
        Cantera::showErrors(std::cout);
        return 1;
    }
    return 0;
}
