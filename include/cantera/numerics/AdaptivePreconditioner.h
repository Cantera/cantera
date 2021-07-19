/**
 *  @file AdaptivePreconditioner.h
 *   Declarations for the class AdaptivePreconditioner
 *   which is a child class of PreconditionerBase for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "float.h"

//! Flag to indicate adaptive preconditioner is set
const int ADAPTIVE_MECHANISM_PRECONDITIONER = 1;

namespace Cantera //Making ASP apart of Cantera namespace
{

  class AdaptivePreconditioner : public PreconditionerBase
    {
    protected:
        //! @param m_matrix is the container that is the sparse preconditioner
        Eigen::SparseMatrix<double> m_matrix;

        //! @param m_identity is the container that is the sparse preconditioner
        Eigen::SparseMatrix<double> m_identity;

        //! @param m_solver is the solver used in solving the linear
        //! system
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solver;

        //! @param m_nonzeros is a member variable for the reserved nonzero elements in m_matrix
        size_t m_nonzeros; // initialize to 0 to prevent failures

        //! @param m_threshold a double value to selectively fill the matrix structure based on this threshold
        double m_threshold = DBL_EPSILON; // default

        //! @param m_current_start an index value for the starting index
        //! of the current reactor in the state
        size_t m_current_start;

        //! @param m_atol absolute tolerance of the ODE solver
        double m_atol = 0;

        // TODO: double check beta and time step
        //! @param m_dampeningParam a value, beta, between zero and 1 to
        //! damp the preconditioner. By default it is undamped, beta = 1.
        double m_dampeningParam = 1;

        //! @param m_timestep the delta t value used in gamma = (t_n-t_(n-1))*beta
        double m_timestep;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW //Required for mis-alignment of EIGEN matrix
        AdaptivePreconditioner(/* args */){}; //Default constructor
        virtual ~AdaptivePreconditioner(){}; //Default destructor
        AdaptivePreconditioner(const AdaptivePreconditioner &externalPrecon); //Copy constructor

        //! This function determines rate law derivatives of species
        //! with respect to other species specifically it determines the
        //! derivatives of the rate laws of all species with respect to
        //! other species in terms of moles.
        //! @param *reactor A pointer to the current reactor being used
        //! for preconditioning
        void SpeciesSpeciesDerivatives(Reactor* reactor);

        //! This function is a subfunction of SpeciesDerivatives that
        //! gets the species w.r.t species derivatives for each reaction
        //! @param reactor pointer to the current reactor being
        //! preconditioned
        //! @param dependent a pointer to a Composition of the species
        //! rate laws that will be differentiated
        //! @param independent a pointer to a Composition of the species
        //! and the derivatives are taken with respect too
        //! @param rateLawDerivatives a pointer to vector data that
        //! temporarily stores the derivatives
        //! @param concentrations an array of concentrations of the current reactor
        //! @param kDirection is reaction rate constants for either forward or reverse
        inline void updateRateLawDerivatives(Reactor *reactor, Composition *dependent, Composition *independent, double* rateLawDerivatives, double* concentrations, double kDirection);


        //! This function determines derivatives of Species and Temperature with respect to Temperature for jacobian preconditioning with a finite difference.
        //! @param reactor A pointer to the current reactor being precondition
        //! @parama t A double value of the current time
        //! @param y A pointer to the current state passed from CVODES
        //! @param ydot A pointer to the current state derivatives
        //! passed from CVODES
        //! @param params A double pointer to sensitivty parameters.
        virtual void TemperatureDerivatives(IdealGasConstPressureReactor* reactor, double t, double* y, double* ydot, double* params);

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

        //! This function checks if there was an error with eigen and
        //! throws it if so.
        void preconditionerErrorCheck();

        //! This function returns the current type of preconditioner as an integer
        virtual size_t getPreconditionerType(){return ADAPTIVE_MECHANISM_PRECONDITIONER;};

        //! Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //! @param reactors A vector pointer of Reactor pointers in the network
        //! @param output a double pointer to the vector (array) to store inv(A)*b
        //! @param rhs_vector a double pointer to the vector (array) multiplied by inv(A)
        //! @param size a unsigned length of the vectors
        virtual void solve(std::vector<Reactor*>* reactors,std::vector<size_t>* reactorStart,double* output, double *rhs_vector,size_t size);

        //! This function performs the setup of the preconditioner for a Reactor and should be overloaded for each different reactor time
        //! @param reactors A vector pointer of reactor pointers in the network
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void setup(std::vector<Reactor*>* reactors,std::vector<size_t>* reactorStart, double t, double* y, double* ydot, double* params);

        //! This function is called during setup for any processes that need to be completed prior to setup functions used in sundials.
        //! @param dims A pointer to a dimensions array
        //! @param atol absolute tolerance value of ODE solver
        virtual void initialize(std::vector<size_t> *dims, double atol);

        //! This function is called during setup to clean up previous
        //! setup data
        virtual void reset();

        //! Function used to get index start of the current reactor variable
        virtual double getReactorStart();

        //! Function used to get a specific element of the matrix structure
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        virtual double getElement(size_t row, size_t col); //get element

        //! Function used to return compressed version of the matrix structure
        virtual Eigen::SparseMatrix<double>* getMatrix();

        //! Use this function to get the threshold value for setting elements
        virtual double getThreshold();

        //! Use this function to return the used absolute tolerance
        virtual double getAbsoluteTolerance();

        //! Use this function to get the dampening parameter
        virtual double getDampeningParameter();

        //! Use this function to get the current time step to be used in
        //! preconditioning
        virtual double getTimeStep();

        //! Function used to set index start of the current reactor variable
        virtual void setReactorStart(size_t reactorStart);

        //! Function used to set a specific element of the matrix structure
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        //! @param element double value to be inserted into matrix structure
        virtual void setElement(size_t row, size_t col, double element);//set element

        //! Function used to set compressed version of the matrix structure
        //! @param sparseMatrix a SUNMatrix pointer to a type of SUNMatrix
        virtual void setMatrix(Eigen::SparseMatrix<double>* sparseMatrix);

        //! Use this function to set the threshold value to compare elements against
        //! @param threshold double value used in setting by threshold
        virtual void setThreshold(double threshold);

        //! Use this function to set the absolute tolerance in the
        //! solver outside of the network initialization
        //! @param atol the specified tolerance
        virtual void setAbsoluteTolerance(double atol);

        //! Use this function to set the dampening parameter
        //! @param dampeningParam the desired dampening parameter
        //! between zero and one.
        virtual void setDampeningParameter(double dampeningParam);

        //! Use this function to set the timestep
        //! @param timestep to be used in preconditioning
        virtual void setTimeStep(double timestep);

        //!Use this function to transform Jacobian into preconditioner
        virtual void transformJacobianToPreconditioner();

        //! Function used to complete individual reactor setups
        //! @param reactor A IdealGasConstPressureReactor pointer
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void reactorLevelSetup(IdealGasConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params);

        //! @param reactor - the contents of this reactor will be printed
        inline void printReactorComponents(Reactor* reactor);

        //! Print preconditioner contents
        void printPreconditioner();

        //! Overloading of the == operator to compare values inside preconditioners
        //! @param externalPrecon - == comparison with this object
        bool operator== (const AdaptivePreconditioner &externalPrecon);

        //! Overloading of the = operator to copy one preconditioner to another
        //! @param externalPrecon. Preconditioner becoming this object
        void operator= (const AdaptivePreconditioner &externalPrecon);
      };
}
#endif