#include "cantera/numerics/AdaptivePreconditioner.h"
#include "float.h"
#include <iostream>

namespace Cantera
{
    /**
     * 
     * AdaptivePreconditioner implementations
     * 
     * **/

    void AdaptivePreconditioner::setDimensions(size_t nrows,size_t ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
    }

    void AdaptivePreconditioner::setElement(size_t row,size_t col,double element)
    {
        this->matrix.coeffRef(row,col)=element;
    }

    double AdaptivePreconditioner::getElement(size_t row,size_t col)
    {
        warn_user("AdaptivePreconditioner::getElement","getElement not properly implemented yet, returning 1.0");
        return 1.0;//FIXME
    }

    Eigen::SparseMatrix<double>* AdaptivePreconditioner::getMatrix()
    {
        return &(this->matrix);
    }

    void AdaptivePreconditioner::setMatrix(Eigen::SparseMatrix<double> *sparseMat)
    {
        this->matrix=*(sparseMat);
    }

    void AdaptivePreconditioner::solve(std::vector<Reactor*>* reactors, std::vector<size_t>* reactorStart, double* output, double* rhs_vector,size_t size)
    {   
        double *rhs_vector_temp = new double[size];
        //rhs_vector is currently the mass fractions that come in 
        for (size_t n = 0; n < reactors->size(); n++) 
        {
            ForwardConversion(reactors->at(n),rhs_vector_temp,rhs_vector,reactorStart->at(n));
        }
        //Compressing sparse matrix structure
        this->matrix.makeCompressed();
        //Creating vectors in the form of //Ax=b
        Eigen::Map<Eigen::VectorXd> bVector(rhs_vector_temp,size);
        // Eigen::Map<Eigen::VectorXd> xVector(x,size);
        Eigen::VectorXd xVector;
        //Create eigen solver object
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        // Initialize solver object
        solver.compute(this->matrix);
        this->checkEigenError("AdaptivePreconditioner::solve (Eigen Decomposition)",solver.info());
        //Solve for xVector
        xVector = solver.solve(bVector);
        int flag = this->checkEigenError("AdaptivePreconditioner::solve (Eigen Solve)",solver.info());
        //Copy x vector to x
        if(flag==0) //If there are no errors copy x-vector
        {
            Eigen::VectorXd::Map(output,xVector.rows())= xVector;
        }
        //Convert output back
        for (size_t n = 0; n < reactors->size(); n++) 
        {
            BackwardConversion(reactors->at(n),output,reactorStart->at(n));
        }
        delete[] rhs_vector_temp;
    }

    void AdaptivePreconditioner::setup(std::vector<Cantera::Reactor*>* reactors, std::vector<size_t>* reactorStart, double t, double* y, double* ydot, double* params)
    {
        for (size_t n = 0; n < reactors->size(); n++) {
            (reactors->at(n))->acceptPreconditioner(this,reactorStart->at(n),t,y,ydot,params);
        }
    }
    
    void AdaptivePreconditioner::reactorLevelSetup(IdealGasConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params)
    {

        StateMap *reactorStateMap = new StateMap; //Index map of idealGasConstPressureReactor
        
        reactorStateMap->operator[]("start") = reactorStart;
        reactorStateMap->operator[]("numberOfSpecies") = (reactor->getKineticsMgr())->nTotalSpecies();
        reactorStateMap->operator[]("species") = reactorStart+(reactor->neq()-reactorStateMap->operator[]("numberOfSpecies")); //Starting idx for species
        
        for (size_t i = 0; i < reactor->neq(); i++)
        {
            reactorStateMap->operator[](reactor->componentName(i))=i;
        }
        
        double* rateLawDervs = new double[(reactor->getKineticsMgr())->nTotalSpecies()*(reactor->getKineticsMgr())->nTotalSpecies()]; //For later use of rate law derivatives

        SpeciesSpeciesDerivatives(reactor,reactorStateMap,rateLawDervs); //SpeciesDerivatives
        
        TemperatureDerivatives(reactor,reactorStateMap,t,y,ydot,rateLawDervs,params); //Temperature Derivatives
        
        NoPrecondition(reactorStateMap,"mass"); //No precondition on mass variable
        
        delete[] rateLawDervs;
        delete reactorStateMap;
        // std::cout<<Eigen::MatrixXd(this->matrix)<<std::endl;
    }

    void AdaptivePreconditioner::reactorLevelSetup(Reactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params)
    {
        throw CanteraError("AdaptivePrecondtioner:reactorLevelSetup",reactor->typeStr()+" setup not implemented.");
    }

    void AdaptivePreconditioner::reactorLevelSetup(IdealGasReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params)
    {
        throw CanteraError("AdaptivePrecondtioner:reactorLevelSetup",reactor->typeStr()+" setup not implemented.");
    }

    void AdaptivePreconditioner::reactorLevelSetup(ConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params)
    {
        throw CanteraError("AdaptivePrecondtioner:reactorLevelSetup",reactor->typeStr()+" setup not implemented.");
    }

    void AdaptivePreconditioner::initialize(size_t nrows,size_t ncols)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
        this->nonzeros = this->dimensions[0]*this->dimensions[1]/2; //Reserves up to half the total spaces
        this->matrix.resize(nrows,ncols);
        this->matrix.reserve(this->nonzeros);
    }

    void AdaptivePreconditioner::reset()
    {
        //Do any reset stuff here
        this->matrix.setZero(); //Set all elements to zero
        this->matrix.makeCompressed(); //Compress matrix
        this->matrix.reserve(this->nonzeros); //Reserve space potentially needed
    }

    void AdaptivePreconditioner::ForwardConversion(Reactor *currReactor, double *tempState, double *rhs, size_t reactorStart)
    {
        ThermoPhase *thermo = currReactor->getThermoMgr();
        double currMass =  currReactor->mass();
        size_t nStateVars = currReactor->neq()-thermo->nSpecies(); 
        //Transferring unchanged parameters for each reactor
        for (size_t i = 0; i < nStateVars; i++)
        {   
            size_t globalIndex = reactorStart+i;
            tempState[globalIndex] = rhs[globalIndex];
        }
        //Adjusting mass fraction parameters for each reactor to moles for AJP
        double *molecularWeights = new double[thermo->nSpecies()];
        thermo->getMolecularWeights(molecularWeights);
        for (size_t i = 0; i < thermo->nSpecies(); i++)
        {
            size_t globalIndex = reactorStart+i+nStateVars;
            tempState[globalIndex] = currMass*molecularWeights[i]*rhs[globalIndex];     
        }
        delete[] molecularWeights;
    }

    void AdaptivePreconditioner::BackwardConversion(Reactor *currReactor, double *output, size_t reactorStart)
    {
            ThermoPhase *thermo = currReactor->getThermoMgr();
            double currMass =  currReactor->mass();
            size_t nStateVars = currReactor->neq()-thermo->nSpecies(); 
            //Do nothing to unchanged parameters
            //Convert moles back to mass fractions
            double *molecularWeights = new double[thermo->nSpecies()];
            thermo->getMolecularWeights(molecularWeights);
            for (size_t i = 0; i < thermo->nSpecies(); i++)
            {
                size_t globalIndex = reactorStart+i+nStateVars;
                output[globalIndex] *= 1/(molecularWeights[i]*currMass);     
            }
            delete[] molecularWeights;
    }

    void AdaptivePreconditioner::SpeciesSpeciesDerivatives(Reactor* reactor,StateMap* stateMap,double* rateLawDerivatives)
    {
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        //Getting thermophase object for access to concentrations and species data
        ThermoPhase* thermo=reactor->getThermoMgr();
        //Compositions for reactants and products
        Composition reactants;
        Composition products;
        //Important sizes to the determination of values
        size_t numberOfReactions = kinetics->nReactions();
        size_t numberOfSpecies = stateMap->operator[]("numberOfSpecies");
        size_t speciesStart  = stateMap->operator[]("species"); //Starting idx for species 
        //Array pointers for data that is reused
        double* kForward = new double[numberOfReactions];
        double* kBackward = new double[numberOfReactions];
        //Concentrations of species
        double* concentrations = new double[numberOfSpecies];
        //Getting species concentrations
        thermo->getConcentrations(concentrations);
        //Getting forward rate constants for calcs
        kinetics->getFwdRateConstants(kForward); 
        //Getting reverse rate constants for calcs
        kinetics->getRevRateConstants(kBackward); 
        //Getting forward and reverse rates of progress
        for (size_t r = 0; r < numberOfReactions; r++)
        {
            std::shared_ptr<Reaction> currentReaction=kinetics->getReactionPtr(r); //shared_ptr for current reaction in finding Jacobian
            //Loop through reactants in current reaction
            reactants = currentReaction->reactants;
            products = currentReaction->products;
            this->GetRateOfProgress(reactants,stateMap,rateLawDerivatives,concentrations,kForward[r],reactor->volume(),numberOfSpecies);
            this->GetRateOfProgress(products,stateMap,rateLawDerivatives,concentrations,-1*kBackward[r],reactor->volume(),numberOfSpecies); //Multiply by negative one to change direction
        }
        
        //Adding to preconditioner
        //d(w)/dn_j
        for (size_t j = 0; j < numberOfSpecies; j++) // column
            {
            for (size_t i = 0; i < numberOfSpecies; i++) //row
            {  
            size_t idx = j+i*numberOfSpecies; //Getting flattened index
            this->setElementByThreshold(i+speciesStart,j+speciesStart,reactor->volume()*rateLawDerivatives[idx]);//Add by threshold 
            }
        }
        //Deleting appropriate pointers
        delete[] kForward;
        delete[] kBackward;
        delete[] concentrations;
    }

    inline void AdaptivePreconditioner::GetRateOfProgress(std::map<std::string, double> comp, StateMap* stateMap, double* omega, double* concentrations, double k_direction, double volume, size_t numberOfSpecies)
    { 
        //flattened index for derivatives
        size_t oidx; //index for omega
        size_t sidx; //index for species
        size_t speciesStart = stateMap->operator[]("species"); //Adjustment based on reactor, network, and stateMap to get Omega index
        double dRdn; //temporary value for rate derivative
        for (std::map<std::string,double>::iterator iter1 = comp.begin(); iter1 != comp.end(); iter1++) //Independent variable -- column
        {
            for (std::map<std::string,double>::iterator iter2 = comp.begin(); iter2 != comp.end(); iter2++) //Dependent variable -- row
            {
            //Get index for current omega
            oidx = stateMap->operator[](iter1->first)-speciesStart+(stateMap->operator[](iter2->first)-speciesStart)*numberOfSpecies;
            dRdn=1; //Set dRdn to one so multiplication isn't zero unless a coefficient is zero
            // Loop through species to get rate derivative
            for (std::map<std::string,double>::iterator iter3 = comp.begin(); iter3 != comp.end(); iter3++) //Dependent variable
            {
                sidx = stateMap->operator[](iter3->first)-speciesStart;
                if (iter3->first == iter1->first)
                {
                dRdn *= std::pow(iter3->second*concentrations[sidx],iter3->second-1); //derivative
                }
                else
                {
                dRdn *= std::pow(concentrations[sidx],iter3->second); //Not derivative
                }
            }
            omega[oidx] += k_direction*dRdn/volume; //Updating omega derivative as is necessary
            }
        }
    } 

    void AdaptivePreconditioner::TemperatureDerivatives(IdealGasConstPressureReactor* reactor,StateMap* stateMap, double t, double* y, double* ydot, double* rateLawDerivatives, double* params)
    {   
        //Getting kinetics object for access to reactions
        Kinetics* kinetics=reactor->getKineticsMgr();
        ThermoPhase* thermo=reactor->getThermoMgr();
        //Important sizes to the determination of values
        size_t numberOfSpecies = kinetics->nTotalSpecies();
        size_t speciesStart  = stateMap->operator[]("species"); //Starting idx for species 
        size_t tempIndex = stateMap->operator[]("temperature")+stateMap->operator[]("start");
        
        //Array pointers for data that is reused
        //net production rates (omega dot)
        double* netProductionRatesNext = new double[numberOfSpecies];
        double* netProductionRatesCurrent = new double[numberOfSpecies];
        
        //Getting perturbed state
        //Perturbation for finite difference of temperature
        double deltaTemp = y[tempIndex]*(std::sqrt(DBL_EPSILON));
        thermo->setTemperature(y[tempIndex]+deltaTemp);
        kinetics->getNetProductionRates(netProductionRatesNext);
        double TDotNext = reactor->evaluateEnergyEquation(t,y,ydot,params)/(thermo->cp_mass()*reactor->mass()); //Perturbed internal energy
        
        //Getting current state
        thermo->setTemperature(y[tempIndex]); //Setting temperature back to correct value
        kinetics->getNetProductionRates(netProductionRatesCurrent);
        double TDotCurrent = reactor->evaluateEnergyEquation(t,y,ydot,params)/(thermo->cp_mass()*reactor->mass()); //Current internal energy

        /**
         * Temp Rate Derivatives w.r.t Temp
         * d T_dot/dT
         **/
        this->setElementByThreshold(tempIndex,tempIndex,(TDotNext-TDotCurrent)/deltaTemp);
        /**
         * Production Rate Derivatives w.r.t Temp
         * d omega_dot_j/dT
         **/
        //convert kmol/m^3/s to kmol/s by multiplying volume and deltaTemp
        for (size_t j = 0; j < numberOfSpecies; j++) //column
        {   
            this->setElementByThreshold(j+speciesStart,tempIndex,(netProductionRatesNext[j]-netProductionRatesCurrent[j])/deltaTemp); //Add by threshold specTempDerivative
        }
        //Deleting appropriate pointers
        delete[] netProductionRatesNext;
        delete[] netProductionRatesCurrent;

        /**
         * Temp Rate Derivatives w.r.t Species
         * d T_dot/dnj
         **/
        double* enthalpy = new double[numberOfSpecies];
        double* specificHeat = new double[numberOfSpecies];
        double* netProductionRates = new double[numberOfSpecies];
        double* concentrations = new double[numberOfSpecies];
        //Getting species concentrations
        thermo->getConcentrations(concentrations);
        thermo->getPartialMolarEnthalpies(enthalpy);
        thermo->getPartialMolarCp(specificHeat);
        kinetics->getNetProductionRates(netProductionRates);
        //Getting perturbed changes w.r.t temperature
        for (size_t j = 0; j < numberOfSpecies; j++) //Spans columns
        {
            double CkCpkSum = 0;
            double hkwkSum = 0;
            double hkdwkdnjSum = 0;
            for (size_t k = 0; k < numberOfSpecies; k++) //Spans rows
            {   
                int idx = j+k*numberOfSpecies; //Getting flattened index - j to remain same and k to change. This means moving down a row.
                hkwkSum += enthalpy[k]*netProductionRates[k];
                hkdwkdnjSum += enthalpy[k]*rateLawDerivatives[idx];
                CkCpkSum += concentrations[k]*specificHeat[k];
            }
            //Set appropriate colume of preconditioner
            this->setElementByThreshold(tempIndex,j+speciesStart,(-hkdwkdnjSum*CkCpkSum+specificHeat[j]/reactor->volume()*hkwkSum)/(CkCpkSum*CkCpkSum));
        }
        delete[] enthalpy;
        delete[] specificHeat;
        delete[] netProductionRates;
        delete[] concentrations;


    }


    int AdaptivePreconditioner::checkEigenError(std::string method, size_t info)
    {   
        int flag = 0;
        if(info!=Eigen::Success) 
        {   
            
            std::string error="Failure: ";
            if(info==Eigen::NumericalIssue)
            {
                error+="NumericalIssues";
                flag+=1;
            }
            else if(info==Eigen::NoConvergence)
            {
                error+="NoConvergence";
                flag+=2;
            }
            else if(info==Eigen::InvalidInput)
            {
                error+="InvalidInput";
                flag+=3;
            }
            else
            {
                error+="Unknown";
                flag+=4;
            }
            warn_user(method,error);
            // throw CanteraError(method,error); //Error throw causing segfault?
        }
        return flag;
    }

    void AdaptivePreconditioner::NoPrecondition(StateMap* stateMap, std::string key)
    {   
        size_t idx = stateMap->operator[](key)+stateMap->operator[]("start");
        this->setElement(idx,idx,1); //setting key variable element of preconditioner equal to 1
    }

    /*
        Other functions used in preconditioner functions but not directly related to a state variable
    */
    
    inline void AdaptivePreconditioner::printReactorComponents(Reactor* reactor)
    {
        for (size_t i = 0; i < reactor->neq(); i++)
        {
        std::cout<<reactor->componentName(i)<<std::endl;
        }
    }

}