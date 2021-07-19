/**
 *  @file PreconditionerBase.h
 *   Declarations for the class PreconditionerBase which is a virtual base class for preconditioning systems.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef PRECONDITIONERBASE_H
#define PRECONDITIONERBASE_H

//Cantera Imports //Need reactor includes for reactorLevelSetup
#include "cantera/base/ctexceptions.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/FlowReactor.h"

//! Flag to indicate preconditioner is not set
const int PRECONDITIONER_NOT_SET = 0;

namespace Cantera//Making ASP apart of Cantera namespace
{
  class PreconditionerBase
  {
    protected:
        //! @param dimensions an size_t pointer to store dimensions
        std::vector<size_t> m_dimensions;
    public:
        PreconditionerBase(/* args */){} //default constructor
        ~PreconditionerBase(){} //default destructor
        //! This function returns zero for preconditioner not set
        virtual size_t getPreconditionerType(){return PRECONDITIONER_NOT_SET;};

        //! Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //! @param reactors A vector pointer of reactor pointers in the network
        //! @param output a double pointer to the vector (array) to store inv(A)*b
        //! @param rhs_vector a double pointer to the vector (array) multiplied by inv(A)
        //! @param size a unsigned length of the vectors
        virtual void solve(std::vector<Reactor*>* reactors, std::vector<size_t>* reactorStart, double* output, double *rhs_vector,size_t size){
            throw CanteraError("Cantera::PrecondtionerBase::solve","solve not implemented for PreconditionerBase, please use a subclass.");
        };

        //! This function performs the setup of the preconditioner for the specified reactor type and should be overloaded for each different reactor time
        //! @param reactors A vector to reactor objects in the network
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void setup(std::vector<Reactor*>* reactors, std::vector<size_t>* reactorStart, double t, double* y, double* ydot, double* params){
            throw CanteraError("Cantera::PrecondtionerBase::setup","setup not implemented for PreconditionerBase, please use a subclass.");
        };

        //! This function is called during setup for any processes that need to be completed prior to setup functions
        //! @param dims A pointer to a dimensions array
        //! @param atol absolute tolerance of the ODE solver
        virtual void initialize(std::vector<size_t> *dims, double atol){
            throw CanteraError("Cantera::PrecondtionerBase::initialize","initialize not implemented for PreconditionerBase, please use a subclass.");
        };

        //! This function is called during setup for any processes that need to be completed post to setup functions
        virtual void reset(){
            throw CanteraError("Cantera::PrecondtionerBase::reset","reset not implemented for PreconditionerBase, please use a subclass.");
        };

        //! Function used to set a specific element of the matrix structure
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        //! @param element double value to be inserted into matrix structure
        virtual void setElement(size_t row, size_t col, double element){
            throw CanteraError("Cantera::PrecondtionerBase(setElement","setElement is not defined for PreconditionerBase, please use a subclass.");
        };

        //! Function used to get a specific element of the matrix structure
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        virtual double getElement(size_t row, size_t col){
            throw CanteraError("Cantera::PrecondtionerBase(getElement","getElement is not defined for PreconditionerBase, please use a subclass.");
        };

        //! Function used to set a specific element of the matrix structure
        //! @param timestep double value to be used when solving for the
        //! preconditioner.
        virtual void setTimeStep(double timestep){
            throw CanteraError("Cantera::PrecondtionerBase(setTimeStep","setTimeStep is not defined for PreconditionerBase, please use a subclass.");
        };

        //! Function used to get a timestep value
        virtual double getTimeStep(){
            throw CanteraError("Cantera::PrecondtionerBase(getTimeStep","getTimeStep is not defined for PreconditionerBase, please use a subclass.");
        };

        //! Function used to complete individual reactor setups
        //! @param reactor A Reactor pointer
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void reactorLevelSetup(Reactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params){
        throw CanteraError("Cantera::PrecondtionerBase::reactorLevelSetup","reactorLevelSetup not implemented for specified reactor or preconditioner type.");
        };

        //! Function used to complete individual reactor setups
        //! @param reactor A IdealGasConstPressureReactor pointer
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void reactorLevelSetup(IdealGasConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params){
        throw CanteraError("Cantera::PrecondtionerBase::reactorLevelSetup","Cantera::PrecondtionerBase::reactorLevelSetup","reactorLevelSetup not implemented for specified reactor or preconditioner type.");
        };

        //! Function used to set dimensions of the preconditioner
        //! @param dims A pointer to an array of the dimensions
        void setDimensions(std::vector<size_t> *dims)
        {
            this->m_dimensions.clear();
            for (auto it = dims->begin(); it != dims->end(); ++it)
            {
                this->m_dimensions.push_back(*it);
            }
        }

        //! Function to return pointer to dimensions
        std::vector<size_t>* getDimensions()
        {
            return &(this->m_dimensions);
        }
  };
}
#endif
