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

//Preconditioner type
const int PRECONDITIONER_NOT_SET = 0;

namespace Cantera//Making ASP apart of Cantera namespace
{
  class PreconditionerBase
  {
    private:

    protected:
        //@param threshold a double value to selectively fill the matrix structure based on this threshold
        double threshold=10e-16; //default
        //@param dimensions an size_t pointer to store dimensions
        size_t dimensions[2];
    public:
        PreconditionerBase(/* args */){}
        PreconditionerBase(const PreconditionerBase &precBase){*this=precBase;} //Copy constructor
        virtual ~PreconditionerBase(){} //destructor
        virtual size_t getPreconditionerType(){return PRECONDITIONER_NOT_SET;};
        /*

          Reactor Setup & Solve Functions

        */
       //!Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //!@param reactors A vector pointer of reactor pointers in the network
        //@param output a double pointer to the vector (array) to store inv(A)*b
        //@param rhs_vector a double pointer to the vector (array) multiplied by inv(A)
        //@param size a unsigned length of the vectors
        virtual void solve(std::vector<Reactor*>* reactors, std::vector<size_t>* reactorStart, double* output, double *rhs_vector,size_t size){
          throw CanteraError("Cantera::PrecondtionerBase::solve","solve not implemented for PreconditionerBase, please use a subclass.");
        };

        //! This function performs the setup of the preconditioner for the specified reactor type and should be overloaded for each different reactor time
        //!@param reactors A vector to reactor objects in the network
        //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void setup(std::vector<Reactor*>* reactors, std::vector<size_t>* reactorStart, double t, double* y, double* ydot, double* params){
          throw CanteraError("Cantera::PrecondtionerBase::setup","setup not implemented for PreconditionerBase, please use a subclass.");
        };

        //!This function is called during setup for any processes that need to be completed prior to setup functions
        //! e.g. dynamic memory allocation
        virtual void initialize(size_t nrows,size_t ncols){
          throw CanteraError("Cantera::PrecondtionerBase::initialize","initialize not implemented for PreconditionerBase, please use a subclass.");
        };

        //!This function is called during setup for any processes that need to be completed post to setup functions
        //! e.g. dynamic memory allocation
        virtual void reset(){
          throw CanteraError("Cantera::PrecondtionerBase::reset","reset not implemented for PreconditionerBase, please use a subclass.");
        };

        //!Function used to set a specific element of the matrix structure
        //!@param row size_t specifying the row location
        //!@param col size_t specifying the column location
        //!@param element double value to be inserted into matrix structure
        virtual void setElement(size_t row, size_t col, double element)
        {
          throw CanteraError("Cantera::PrecondtionerBase::setElement","setElement not implemented for PreconditionerBase, please use a subclass.");
        }; //set element

        //!Function used to get a specific element of the matrix structure
        //!@param row size_t specifying the row location
        //!@param col size_t specifying the column location
        virtual double getElement(size_t row, size_t col){
          throw CanteraError("Cantera::PrecondtionerBase::getElement","getElement not implemented for PreconditionerBase, please use a subclass.");
        }; //get element

        //Other preconditioner functions
        //!Use this function to get the threshold value for setting elements
        virtual double getThreshold();

        //!Use this function to set the threshold value to compare elements against
        //!@param threshold double value used in setting by threshold
        virtual void setThreshold(double threshold);

        //!Use this function to set an element by the threshold
        //!@param row size_t specifying the row location
        //!@param col size_t specifying the column location
        //!@param element double value to be inserted into matrix structure
        virtual void setElementByThreshold(size_t row,size_t col, double element);

        //!Function used to set the dimensions of and construct the matrix structure - required for initialization and use of the class
        //!@param nrows size_t number of rows in the structure
        //!@param ncols size_t nubmer of columns in the structure
        //!@param otherData void* for passing other data necessary for subclasses to initialize the matrix structure
        virtual void setDimensions(size_t nrows,size_t ncols);

        //!Function to return the dimensions of the matrix structure
        virtual size_t* getDimensions();

        /**
         *
         * Reactor Level Functions
         *
         **/

        //!Function used to complete individual reactor setups
          //!@param reactor A Reactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(Reactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params){
          throw CanteraError("Cantera::PrecondtionerBase::reactorLevelSetup","reactorLevelSetup not implemented for PreconditionerBase, please use a subclass.");
        };

          //!Function used to complete individual reactor setups
          //!@param reactor A IdealGasConstPressureReactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(IdealGasConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params){
          throw CanteraError("Cantera::PrecondtionerBase::reactorLevelSetup","reactorLevelSetup not implemented for PreconditionerBase, please use a subclass.");
        };

          //!Function used to complete individual reactor setups
          //!@param reactor A IdealGasReactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(IdealGasReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params){
          throw CanteraError("Cantera::PrecondtionerBase::reactorLevelSetup","reactorLevelSetup not implemented for PreconditionerBase, please use a subclass.");
        };

          //!Function used to complete individual reactor setups
          //!@param reactor A IdealGasReactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(ConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params){
          throw CanteraError("Cantera::PrecondtionerBase::reactorLevelSetup","reactorLevelSetup not implemented for PreconditionerBase, please use a subclass.");
        };

  };
}
#endif
