#include "cantera/numerics/PreconditionerBase.h"

namespace Cantera
{
    /**
     * 
     * Preconditioner Implementation
     * 
     * **/
    double PreconditionerBase::getThreshold()
    {
        return this->threshold;
    }

    void PreconditionerBase::setThreshold(double threshold)
    {
        this->threshold = threshold;
    }

    void PreconditionerBase::setElementByThreshold(size_t row,size_t col, double element)
    {
        if (this->threshold<element)
        {
            this->setElement(row,col,element);
        }
    }

    void PreconditionerBase::setDimensions(size_t nrows,size_t ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
    }

    size_t* PreconditionerBase::getDimensions()
    {
        return &(this->dimensions[0]);
    }
}
