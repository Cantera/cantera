//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/AdaptivePreconditioner.h"
#include "float.h"
#include <iostream>

namespace Cantera
{

    AdaptivePreconditioner::AdaptivePreconditioner(const AdaptivePreconditioner &externalPrecon)
    {
        m_threshold = externalPrecon.m_threshold;
        m_matrix = externalPrecon.m_matrix;
        m_nonzeros = externalPrecon.m_nonzeros;
        m_dimensions.clear();
        m_dimensions.push_back(externalPrecon.m_dimensions.at(0));
        m_dimensions.push_back(externalPrecon.m_dimensions.at(1));
    }

    bool AdaptivePreconditioner::operator== (const AdaptivePreconditioner &externalPrecon)
    {
        // only compares internal matrix
        return m_matrix.isApprox(externalPrecon.m_matrix);
    }

    void AdaptivePreconditioner::operator= (const AdaptivePreconditioner &externalPrecon)
    {
        m_threshold = externalPrecon.m_threshold;
        m_matrix = externalPrecon.m_matrix;
        m_nonzeros = externalPrecon.m_nonzeros;
        m_atol = externalPrecon.m_atol;
        m_dimensions.clear();
        m_dimensions.push_back(externalPrecon.m_dimensions.at(0));
        m_dimensions.push_back(externalPrecon.m_dimensions.at(1));
        return;
    }

    void AdaptivePreconditioner::setElement(size_t row, size_t col, double element)
    {
        if (std::abs(element) >= m_threshold)
        {
            m_matrix.coeffRef(row,col) = element;
        }
        else if( row == col)
        {
            m_matrix.coeffRef(row,col) = element;
        }
    }

    double AdaptivePreconditioner::getElement(size_t row, size_t col)
    {
        return m_matrix.coeffRef(row, col);
    }

    Eigen::SparseMatrix<double>* AdaptivePreconditioner::getMatrix()
    {
        return &(m_matrix);
    }

    void AdaptivePreconditioner::setMatrix(Eigen::SparseMatrix<double> *sparseMat)
    {
        m_matrix =* (sparseMat);
    }

    double AdaptivePreconditioner::getThreshold()
    {
        return m_threshold;
    }

    void AdaptivePreconditioner::setThreshold(double threshold)
    {
        m_threshold = threshold;
    }

    void AdaptivePreconditioner::setReactorStart(size_t reactorStart)
    {
        m_current_start = reactorStart;
    }

    double AdaptivePreconditioner::getReactorStart()
    {
        return m_current_start;
    }

    void AdaptivePreconditioner::setAbsoluteTolerance(double atol)
    {
        m_atol = atol;
    }

    double AdaptivePreconditioner::getAbsoluteTolerance()
    {
        return m_atol;
    }

    void AdaptivePreconditioner::setDampeningParameter(double dampeningParam)
    {
        m_dampeningParam = dampeningParam;
    }

    double AdaptivePreconditioner::getDampeningParameter()
    {
        return m_dampeningParam;
    }

    void AdaptivePreconditioner::setTimeStep(double timestep)
    {
        m_timestep = timestep;
    }

    double AdaptivePreconditioner::getTimeStep()
    {
        return m_timestep;
    }

    void AdaptivePreconditioner::solve(std::vector<Reactor*>* reactors, std::vector<size_t>* reactorStart, double* output, double* rhs_vector, size_t size)
    {
        std::vector<double> rhs_vector_temp (size);
        // rhs_vector is currently the mass fractions that come in
        for (size_t n = 0; n < reactors->size(); n++)
        {
            ForwardConversion(reactors->at(n), rhs_vector_temp.data(), rhs_vector, reactorStart->at(n));
        }
        // Compressing sparse matrix structure
        m_matrix.makeCompressed();
        // Creating vectors in the form of Ax=b
        Eigen::Map<Eigen::VectorXd> bVector(rhs_vector_temp.data(), size);
        Eigen::VectorXd xVector;
        // Solve for xVector
        xVector = m_solver.solve(bVector);
        preconditionerErrorCheck();
        // Copy x vector to x
        Eigen::VectorXd::Map(output, xVector.rows()) = xVector;
        // Convert output back
        for (size_t n = 0; n < reactors->size(); n++)
        {
            BackwardConversion(reactors->at(n), output, reactorStart->at(n));
        }
    }

    void AdaptivePreconditioner::setup(std::vector<Cantera::Reactor*>* reactors, std::vector<size_t>* reactorStart, double t, double* y, double* ydot, double* params)
    {
        for (size_t n = 0; n < reactors->size(); n++)
        {
            // Reactor start is not added to y and ydot because the
            // preconditioner is not broken into units like reactors are
            (reactors->at(n))->acceptPreconditioner(this, reactorStart->at(n), t, y, ydot, params);
        }
        // Make into preconditioner as P = (I - gamma * J_bar)
        transformJacobianToPreconditioner();
        // Analyze and factorize
        m_solver.analyzePattern(m_matrix);
        m_solver.factorize(m_matrix);
        preconditionerErrorCheck();
    }

    void AdaptivePreconditioner::reactorLevelSetup(IdealGasConstPressureReactor* reactor, size_t reactorStart, double t, double* y, double* ydot, double* params)
    {
        // Create strictly positive composition
        std::vector<double> yCopy(reactor->neq());
        for (size_t i = 0; i < reactor->neq(); i++)
        {
            yCopy[i] = std::max(y[reactorStart+i], m_atol);
        }
        // Set state to strictly positive composition
        reactor->updateState(yCopy.data());
        // Set preconditioner reactor start to passed value
        setReactorStart(reactorStart);
        // SpeciesDerivatives
        SpeciesSpeciesDerivatives(reactor);
        // Temperature Derivatives
        TemperatureDerivatives(reactor, t, yCopy.data(), ydot, params);
        // Reset state
        reactor->updateState(y+reactorStart);
    }

    void AdaptivePreconditioner::initialize(std::vector<size_t> *dims, double atol)
    {
        m_dimensions.clear();
        for (auto dimIter = dims->begin(); dimIter != dims->end(); ++dimIter)
        {
            m_dimensions.push_back(*dimIter);
        }
        if (m_dimensions[0] != m_dimensions[1])
        {
            throw CanteraError("initialize", "specified matrix dimensions are not square");
        }
        m_nonzeros = m_dimensions[0]; // Reserves the main diagonal spaces
        m_matrix.resize(m_dimensions[0], m_dimensions[1]);
        m_matrix.reserve(m_nonzeros);
        m_atol = atol;
        // Creating sparse identity matrix
        m_identity.resize(m_dimensions[0], m_dimensions[1]);
        m_identity.setIdentity();
    }

    void AdaptivePreconditioner::reset()
    {
        m_matrix.setZero(); // Set all elements to zero
        m_matrix.makeCompressed(); // Compress matrix
        m_matrix.reserve(m_nonzeros); // Reserve space potentially needed
    }

    void AdaptivePreconditioner::ForwardConversion(Reactor *currReactor, double *tempState, double *rhs, size_t reactorStart)
    {
        ThermoPhase *thermo = currReactor->getThermoMgr();
        double currMass =  currReactor->mass();
        size_t nStateVars = currReactor->neq()-thermo->nSpecies();
        // Transferring unchanged parameters for each reactor
        for (size_t i = 0; i < nStateVars; i++)
        {
            size_t globalIndex = reactorStart+i;
            tempState[globalIndex] = rhs[globalIndex];
        }
        // Adjusting mass fraction parameters for each reactor to moles for AJP
        std::vector<double> molecularWeights (thermo->nSpecies());
        thermo->getMolecularWeights(molecularWeights.data());
        for (size_t i = 0; i < thermo->nSpecies(); i++)
        {
            size_t globalIndex = reactorStart + i + nStateVars;
            tempState[globalIndex] = currMass * molecularWeights[i] * rhs[globalIndex];
        }
    }

    void AdaptivePreconditioner::BackwardConversion(Reactor *currReactor, double *output, size_t reactorStart)
    {
        ThermoPhase *thermo = currReactor->getThermoMgr();
        double currMass =  currReactor->mass();
        size_t nStateVars = currReactor->neq() - thermo->nSpecies();
        // Do nothing to unchanged parameters
        // Convert moles back to mass fractions
        std::vector<double> molecularWeights (thermo->nSpecies());
        thermo->getMolecularWeights(molecularWeights.data());
        for (size_t i = 0; i < thermo->nSpecies(); i++)
        {
            size_t globalIndex = reactorStart + i + nStateVars;
            output[globalIndex] *= 1/(molecularWeights[i] * currMass);
        }
    }

    void AdaptivePreconditioner::SpeciesSpeciesDerivatives(Reactor* reactor)
    {
        auto kinetics = reactor->getKineticsMgr();
        auto thermo = reactor->getThermoMgr();
        // Important sizes to the determination of values
        size_t numberOfReactions = kinetics->nReactions();
        size_t numberOfSpecies = kinetics->nTotalSpecies();
        size_t speciesStart = getReactorStart() + thermo->stateSize() - numberOfSpecies;
        // Vectors for data that is reused
        std::vector<double> rateLawDerivatives (numberOfSpecies * numberOfSpecies);
        std::vector<double> kForward (numberOfReactions);
        std::vector<double> kBackward (numberOfReactions);
        std::vector<double> concentrations (numberOfSpecies);
        thermo->getConcentrations(concentrations.data());
        kinetics->getFwdRateConstants(kForward.data());
        kinetics->getRevRateConstants(kBackward.data());
        double volume = reactor->volume();

        for (size_t r = 0; r < numberOfReactions; r++)
        {
            // shared_ptr for current reaction in finding Jacobian
            std::shared_ptr<Reaction> currentReaction = kinetics->getReaction(r);
            Composition reactants = currentReaction->reactants;
            Composition products = currentReaction->products;
            updateRateLawDerivatives(reactor, &reactants, &reactants, rateLawDerivatives.data(), concentrations.data(), -kForward[r]);
            updateRateLawDerivatives(reactor, &products, &reactants, rateLawDerivatives.data(), concentrations.data(), kForward[r]);
            // Calculate other direction if reaction is reversible
            updateRateLawDerivatives(reactor, &reactants, &products, rateLawDerivatives.data(), concentrations.data(), kBackward[r]);
            updateRateLawDerivatives(reactor, &products, &products, rateLawDerivatives.data(), concentrations.data(), -kBackward[r]);
        }
        // Adding to preconditioner
        // d(w)/dn_j
        for (size_t j = 0; j < numberOfSpecies; j++) //  row
            {
            for (size_t i = 0; i < numberOfSpecies; i++) // col
            {
                size_t idx = j + i * numberOfSpecies; // Getting flattened index
                setElement(j + speciesStart, i + speciesStart, volume * rateLawDerivatives[idx]); // Add by threshold
            }
        }
    }

    inline void AdaptivePreconditioner::updateRateLawDerivatives(Reactor *reactor, Composition *dependent, Composition *independent, double* rateLawDerivatives, double* concentrations, double kDirection)
    {
        size_t numberOfSpecies = reactor->getThermoMgr()->nSpecies();
        size_t speciesStart = reactor->getThermoMgr()->stateSize() - numberOfSpecies; // Starting idx for species
        double volume = reactor->volume();
        for (std::map<std::string,double>::iterator iterCol = dependent->begin(); iterCol != dependent->end(); iterCol++) // Dependent variable -- column
        {
            for (std::map<std::string,double>::iterator iterRow = independent->begin(); iterRow != independent->end(); iterRow++) // Independent variable -- row
            {
                // Get flattened index for current omega
                size_t omegaIdx = reactor->componentIndex(iterCol->first) - speciesStart + (reactor->componentIndex(iterRow->first)-speciesStart) * (numberOfSpecies);
                // Current derivative
                double dRdn = kDirection * iterCol->second * iterRow->second / volume * std::pow(concentrations[reactor->componentIndex(iterRow->first)-speciesStart], iterRow->second-1);
                for (std::map<std::string,double>::iterator iterLoop = independent->begin(); iterLoop != independent->end(); iterLoop++) // Dependent variable
                {
                    if (iterRow->first!=iterLoop->first) // Non-derivative terms
                    {
                        dRdn *= std::pow(concentrations[reactor->componentIndex(iterLoop->first)-speciesStart], iterLoop->second);
                    }
                }
                rateLawDerivatives[omegaIdx] += dRdn; // Updating omega derivative
            }
        }
    }

    void AdaptivePreconditioner::TemperatureDerivatives(IdealGasConstPressureReactor* reactor, double t, double* y, double* ydot, double* params)
    {
        auto kinetics = reactor->getKineticsMgr();
        auto thermo = reactor->getThermoMgr();
        // Important sizes to the determination of values
        size_t numberOfSpecies = kinetics->nTotalSpecies();
        // Size of state for current reactor
        size_t stateSize = thermo->stateSize();
        // Starting idx for species in reactor
        size_t speciesStart  = stateSize - numberOfSpecies;
        // Starting idx for species in reactor
        size_t reactorStart = getReactorStart();
        // Temperature Index in reactor
        size_t tempIndex = reactor->componentIndex("temperature");
        // Molecular weights for conversion
        std::vector<double> molecularWeights (thermo->nSpecies());
        thermo->getMolecularWeights(molecularWeights.data());
        // Getting perturbed state for finite difference
        double deltaTemp = y[tempIndex] * (std::sqrt(DBL_EPSILON));
        double reactorMass = reactor->mass();
        // net production rates and enthalpies for each state
        std::vector<double> yNext (stateSize);
        std::vector<double> yDotNext (stateSize);
        std::vector<double> yCurrent (stateSize);
        std::vector<double> yDotCurrent (stateSize);
        // Copy y to current and next
        for (size_t i = 0; i < stateSize; i++)
        {
            yCurrent[i] = y[i];
            yNext[i] = y[i];
        }
        // perturb temperature
        yNext[tempIndex] += deltaTemp;
        // Getting perturbed state
        reactor->updateState(yNext.data());
        reactor->evalEqs(t, yNext.data(), yDotNext.data(), params);
        // Reset and get original state
        reactor->updateState(yCurrent.data());
        reactor->evalEqs(t, yCurrent.data(), yDotCurrent.data(), params);
        // d T_dot/dT
        setElement(tempIndex + reactorStart, tempIndex + reactorStart, reactor->volume() * (yDotNext[tempIndex] - yDotCurrent[tempIndex]) / deltaTemp);
        // d omega_dot_j/dT
        for (size_t j = speciesStart; j < stateSize; j++)
        {
            // Convert dy_j/dt to dN_j/dt by multiply by mass and
            // dividing by MW
            setElement(j+reactorStart, tempIndex+reactorStart, (yDotNext[j] - yDotCurrent[j]) * reactorMass / molecularWeights[j-speciesStart] / deltaTemp);
        }
        // d T_dot/dnj
        std::vector<double> specificHeat (numberOfSpecies);
        std::vector<double> netProductionRates (numberOfSpecies);
        std::vector<double> concentrations (numberOfSpecies);
        std::vector<double> enthalpy (numberOfSpecies);
        // Getting species concentrations
        thermo->getConcentrations(concentrations.data());
        thermo->getPartialMolarCp(specificHeat.data());
        thermo->getPartialMolarEnthalpies(enthalpy.data());
        kinetics->getNetProductionRates(netProductionRates.data());
        // Getting perturbed changes w.r.t temperature
        for (size_t j = 0; j < numberOfSpecies; j++) // Spans columns
        {
            double CkCpkSum = 0;
            double hkwkSum = 0;
            double hkdwkdnjSum = 0;
            for (size_t k = 0; k < numberOfSpecies; k++) // Spans rows
            {
                hkwkSum += enthalpy[k] * netProductionRates[k];
                hkdwkdnjSum += enthalpy[k] * getElement(k+speciesStart, j+speciesStart);
                CkCpkSum += concentrations[k] * specificHeat[k];
            }
            // Set appropriate colume of preconditioner
            setElement(tempIndex + reactorStart, j + speciesStart + reactorStart, (-hkdwkdnjSum * CkCpkSum + specificHeat[j] / reactor->volume() * hkwkSum) / (CkCpkSum * CkCpkSum));
        }
    }

    void AdaptivePreconditioner::transformJacobianToPreconditioner()
    {
        m_matrix = m_identity - m_timestep * m_dampeningParam * m_matrix;
    }

    void AdaptivePreconditioner::printPreconditioner()
    {
        std::cout<<Eigen::MatrixXd(m_matrix)<<std::endl;
    }

    inline void AdaptivePreconditioner::printReactorComponents(Reactor* reactor)
    {
        for (size_t i = 0; i < reactor->neq(); i++)
        {
            std::cout << reactor->componentName(i) << std::endl;
        }
    }

    void AdaptivePreconditioner::preconditionerErrorCheck()
    {
        if (m_solver.info() != Eigen::Success)
        {
            throw CanteraError("AdaptivePreconditioner::solve", m_solver.lastErrorMessage());
        }
    }
}