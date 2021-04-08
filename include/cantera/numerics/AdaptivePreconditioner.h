/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

const int ADAPTIVE_MECHANISM_PRECONDITIONER = 1;

//PreconditionerBase imports
#include "cantera/numerics/PreconditionerBase.h"

//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

//Adaptive preconditioning imports
//Cantera imports
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/Reaction.h"

/**
 * 
 * 
 * Adaptive Mechanism Preconditioning Namespace
 * 
 * 
 * */
namespace Cantera //Making ASP apart of Cantera namespace
{
  //!Typedef used for getting indices based on strings

  typedef std::map<std::string,size_t> StateMap;

  class AdaptivePreconditioner : public PreconditionerBase
    { 
      protected:

        /**
         * 
         * Physics Functions
         * 
         **/
        //! This function determines derivatives of Species with respect to species for jacobian preconditioning;
        //! specifically it determines the derivatives of the rate laws of all species with respect to other species in terms of moles.
        //! @param *preconditioner A pointer to a PreconditionerBase Object for preconditioning the system and storing preconditioner values
        //! @param *reactor A pointer to the current reactor being precondition
        void SpeciesSpeciesDerivatives(Reactor* reactor,StateMap* stateMap,double* rateLawDerivatives);

        //!This function is a subfunction of SpeciesDerivatives that gets the species w.r.t species derivatives for each reaction
        void GetRateOfProgress(std::map<std::string, double> comp, StateMap* stateMap, double* omega, double* concentrations, double k_direction, double volume, size_t numberOfSpecies);
          Eigen::SparseMatrix<double> matrix;
          size_t nonzeros;

        //! This function determines derivatives of Species and Temperature with respect to Temperature for jacobian preconditioning with a finite difference.
        //! @param *preconditioner A pointer to a PreconditionerBase Object for preconditioning the system and storing preconditioner values
        //! @param *reactor A pointer to the current reactor being precondition
        //! @param *ydot A pointer to the current data of ydot passed from CVODES
        //! @param meanSpecificHeat The mean specific heat used based on reactor type
        //! @param index The index location of temperature in the state vector
        virtual void TemperatureDerivatives(IdealGasConstPressureReactor* reactor, StateMap* stateMap, double t, double* y, double* ydot, double* rateLawDerivatives, double* params);
        
        //!This function does not precondition the associated equation by assigning it's preconditioner value to a value of 1
        //!@param row the row index of the variable
        //!@param col the column index of the variable
        void NoPrecondition(StateMap* stateMap, std::string key);

        /**
         * 
         * Other Functions
         * 
         **/
        
        //! This function is used to convert the system from mass fraction to mole fraction for solving the linear system with a mole based jacobian.
        //! @param *reactor A pointer to the current reactor being converted
        //! @param *tempState A double pointer to the temporary state used to solve the linear system
        //! @param *rhs A double pointer provided to preconditioner solve of the initial rhs state
        //! @param reactorStart An size_t for the global index of each species
        void ForwardConversion(Reactor *currReactor, double *tempState, double *rhs, size_t reactorStart);

        //! This function is used to convert the system from moles back to mass fraction after being solved with AMP preconditioner
        //! @param *reactor A pointer to the current reactor being converted
        //! @param *output A double pointer to the output vector of mole values
        //! @param reactorStart An size_t for the global index of each species
        void BackwardConversion(Reactor *currReactor, double *output, size_t reactorStart);

        //! This function determines the rate of progress derivatives given a composition of reactants or products
        int checkEigenError(std::string method, size_t info);

        //! Use this function to print and check reactor components
        inline void printReactorComponents(Reactor* reactor);

      public:
          EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
          AdaptivePreconditioner(/* args */){};
          ~AdaptivePreconditioner(){};
          AdaptivePreconditioner(const AdaptivePreconditioner &preconditioner){*this=preconditioner;} //Copy constructor
          virtual size_t getPreconditionerType(){return ADAPTIVE_MECHANISM_PRECONDITIONER;};

          //!Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
          //!@param reactors A vector pointer of Reactor pointers in the network
          //@param output a double pointer to the vector (array) to store inv(A)*b
          //@param rhs_vector a double pointer to the vector (array) multiplied by inv(A)
          //@param size a unsigned length of the vectors
          virtual void solve(std::vector<Reactor*>* reactors,std::vector<size_t>* reactorStart,double* output, double *rhs_vector,size_t size);

          //! This function performs the setup of the preconditioner for a Reactor and should be overloaded for each different reactor time
          //!@param reactors A vector pointer of reactor pointers in the network
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void setup(std::vector<Reactor*>* reactors,std::vector<size_t>* reactorStart, double t, double* y, double* ydot, double* params);

          //! THis function performs set up 
          //!This function is called during setup for any processes that need to be completed prior to setup functions
          //! e.g. dynamic memory allocation
          virtual void initialize(size_t nrows,size_t ncols);

          //!This function is called during setup for any processes that need to be completed post to setup functions
          //! e.g. dynamic memory allocation
          virtual void reset();

          //!Function used to get a specific element of the matrix structure
          //!@param row size_t specifying the row location
          //!@param col size_t specifying the column location
          virtual double getElement(size_t row, size_t col); //get element

          //!Function used to return compressed version of the matrix structure
          virtual Eigen::SparseMatrix<double>* getMatrix();

          //!Function used to set the dimensions of and construct the matrix structure - required for initialization and use of the class
          //!@param nrows size_t number of rows in the structure
          //!@param ncols size_t nubmer of columns in the structure
          //!@param otherData void* for passing other data necessary for subclasses to initialize the matrix structure
          virtual void setDimensions(size_t nrows,size_t ncols);

          //!Function used to set a specific element of the matrix structure
          //!@param row size_t specifying the row location
          //!@param col size_t specifying the column location
          //!@param element double value to be inserted into matrix structure
          virtual void setElement(size_t row, size_t col, double element);//set element

          //!Function used to set compressed version of the matrix structure
          //!@param sparseMatrix a SUNMatrix pointer to a type of SUNMatrix
          //!@param compress a bool dictating whether or not the set matrix needs compressed or not
          virtual void setMatrix(Eigen::SparseMatrix<double>* sparseMatrix);  

        /**
         * 
         * Reactor Level Functions
         * 
         **/

        //!Function used to complete individual reactor setups
          //!@param reactor A Reactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(Reactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params);

          //!Function used to complete individual reactor setups
          //!@param reactor A IdealGasConstPressureReactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(IdealGasConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params);

          //!Function used to complete individual reactor setups
          //!@param reactor A IdealGasConstPressureReactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(IdealGasReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params);

          //!Function used to complete individual reactor setups
          //!@param reactor A IdealGasConstPressureReactor pointer
          //!@param reactorStart an size_t providing the index location in which the state of the given reactor starts
          virtual void reactorLevelSetup(ConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params);
      };
} 
#endif