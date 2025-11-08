classdef Mixture < handle
    % Mixture Class ::
    %
    %     >> m = ct.Mixture(phases)
    %
    % Class :mat:class:`ct.Mixture` represents mixtures of one or more phases of matter.
    % To construct a mixture, supply a cell array of phases and mole numbers ::
    %
    %     >> gas = ct.Solution('gas.yaml');
    %     >> graphite = ct.Solution('graphite.yaml');
    %     >> mix = ct.Mixture({gas, 1.0; graphite, 0.1});
    %
    % Phases may also be added later using the addPhase method ::
    %
    %     >> water = ct.Solution('water.yaml');
    %     >> mix.addPhase(water, 3.0);
    %
    % Note that the objects representing each phase compute only the
    % intensive state of the phase - they do not store any information
    % on the amount of this phase. Mixture objects, on the other hand,
    % represent the full extensive state.
    %
    % Mixture objects are "lightweight" in the sense that they do not
    % store parameters needed to compute thermodynamic or kinetic
    % properties of the phases. These are contained in the
    % ('heavyweight') phase objects. Multiple mixture objects may be
    % constructed using the same set of phase objects. Each one stores
    % its own state information locally, and synchronizes the phase
    % objects whenever it requires phase properties.
    %
    % :param phases:
    %     Cell array of phases and mole numbers.

    properties (SetAccess = immutable)
        mixID = -1
    end

    properties (SetAccess = public)
        T % Temperature [K].
        P % Pressure [Pa].
    end

    properties (SetAccess = protected)

        phases % Phases in the mixture

        nElements % Number of elements in a mixture.

        nPhases % Number of phases in a mixture.

        nSpecies % Number of species in a mixture.

        chemPotentials % Chemical potentials [J/kmol] of species in the mixture.
    end

    methods
        %% Mixture Class Constructor

        function m = Mixture(phases)
            % Create a :mat:class:`ct.Mixture` object.

            ct.isLoaded(true);

            if nargin > 1
                error('ct.Mixture: wrong number of arguments');
            end

            % Create an empty mixture.
            m.mixID = ct.impl.call('mMix_new');
            m.phases = phases;

            % If phases are supplied, add them
            if nargin == 0
                return
            end

            if ~isa(phases, 'cell')
                error('Enter phases as a cell array.');
            end

            % First column contains the phase objects, and the second
            % column contains the mole numbers of each phase.
            [np, nc] = size(phases);

            % If mole numbers are not defined, default to 1 for all phases.
            if nc == 1
                newColumn = num2cell(ones(np, 1));
                phases = [phases, newColumn];
            elseif nc < 1 || nc > 2
                error('Phases should have one or two columns');
            end

            for n = 1:np
                m.addPhase(phases{n, 1}, phases{n, 2});
            end

            m.T = phases{n, 1}.T;
            m.P = phases{n, 1}.P;

        end

        %% Mixture Class Destructor

        function delete(obj)
            % Delete the :mat:class:`ct.Mixture` object.
            if obj.mixID >= 0
                ct.impl.call('mMix_del', obj.mixID);
            end
        end

        %% Mixture Utility methods

        function display(obj)
            % Display the state of the mixture on the terminal.

            ct.impl.call('mMix_updatePhases', obj.mixID);
            [np, nc] = size(obj.phases);

            for n = 1:np
                s = [sprintf('\n*******************    Phase %d', n) ...
                     sprintf('    ******************************\n\n Moles: %12.6g', ...
                     obj.phaseMoles(n))];
                disp(s);
                display(obj.phases{n, 1});
            end

        end

        function addPhase(obj, phase, moles)
            % Add a phase to a mixture. ::
            %
            %     >> m.addPhase(phase, moles)
            %
            % :param phase:
            %     Instance of class :mat:class:`ct.ThermoPhase` which should be added.
            % :param moles:
            %     Number of moles [kmol] of the ``phase`` to be added to this mixture.

            if ~isa(phase, 'ct.ThermoPhase')
                error('Phase object of wrong type.');
            end

            if ~isa(moles, 'numeric')
                error('Number of moles must be numeric');
            end

            if moles < 0.0
                error('Negative moles');
            end

            ct.impl.call('mMix_addPhase', obj.mixID, phase.tpID, moles);

        end

        %% Mixture Get methods

        function temperature = get.T(obj)
            temperature = ct.impl.call('mMix_temperature', obj.mixID);
        end

        function pressure = get.P(obj)
            pressure = ct.impl.call('mMix_pressure', obj.mixID);
        end

        function n = get.nElements(obj)
            n = ct.impl.call('mMix_nElements', obj.mixID);
        end

        function n = get.nPhases(obj)
            n = ct.impl.call('mMix_nPhases', obj.mixID);
        end

        function n = get.nSpecies(obj)
            n = ct.impl.call('mMix_nSpecies', obj.mixID);
        end

        function mu = get.chemPotentials(obj)
            nsp = obj.nSpecies;
            mu = ct.impl.getArray('mMix_getChemPotentials', nsp, obj.mixID);
        end

        function n = nAtoms(obj, e)
            % Number of atoms of an element in a mixture. ::
            %
            %     >> n = m.nAtoms(e)
            %
            % :param e:
            %     Index of element.
            % :return:
            %     Number of atoms for element e.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.

            n = ct.impl.call('mMix_nPhases', obj.mixID, k - 1, e - 1);
        end

        function n = elementIndex(obj, name)
            % Index of an element. ::
            %
            %     >> n = m.elementIndex(name)
            %
            % :param name:
            %     Name of the element whose index is desired.
            % :return:
            %     Index of element with name ``name``.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.

            n = ct.impl.call('mMix_elementIndex', obj.mixID, name) + 1;
        end

        function n = speciesIndex(obj, k, p)
            % Index of a species in a mixture. ::
            %
            %     >> n = m.speciesIndex(k, p)
            %
            % :param name:
            %     Name of the species whose index is desired.
            % :return:
            %     Index of species with name ``name``.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.

            n = ct.impl.call('mMix_speciesIndex', obj.mixID, k - 1, p - 1) + 1;
        end

        function moles = elementMoles(obj, e)
            % Number of moles [kmol] of an element in a mixture. ::
            %
            %     >> moles = m.elementMoles(e)
            %
            % :param e:
            %    Integer element number.
            % :return:
            %    Moles of element number 'e'. If input 'e' is empty, return
            %    moles of every element in the mixture.

            if nargin == 2
                moles = ct.impl.call('mMix_elementMoles', obj.mixID, e);
            elseif nargin == 1
                nel = obj.nElements;
                moles = zeros(1, nel);

                for i = 1:nel
                    moles(i) = ct.impl.call('mMix_elementMoles', obj.mixID, i-1);
                end

            else
                error('wrong number of arguments');
            end

        end

        function moles = phaseMoles(obj, n)
            % Number of moles [kmol] of a phase in a mixture. ::
            %
            %     >> moles = m.phaseMoles(n)
            %
            % :param n:
            %    Integer phase number.
            % :return:
            %    Moles of element number 'n'. If input 'n' is empty, return
            %    moles of every element in the mixture.

            if nargin == 2
                moles = ct.impl.call('mMix_phaseMoles', obj.mixID, n);
            elseif nargin == 1
                np = obj.nPhases;
                moles = zeros(1, np);

                for i = 1:np
                    moles(i) = ct.impl.call('mMix_phaseMoles', obj.mixID, i-1);
                end

            else
                error('wrong number of arguments');
            end

        end

        function moles = speciesMoles(obj, k)
            % Number of moles [kmol] of a species in a mixture. ::
            %
            %     >> moles = m.speciesMoles(k)
            %
            % :param k:
            %    Integer species number.
            % :return:
            %    Moles of species number 'k'. If input 'k' is empty, return
            %    moles of every species in the mixture

            if nargin == 2
                moles = ct.impl.call('mMix_speciesMoles', obj.mixID, k);
            elseif nargin == 1
                nsp = obj.nSpecies;
                moles = zeros(1, nsp);

                for i = 1:nsp
                    moles(i) = ct.impl.call('mMix_speciesMoles', obj.mixID, i-1);
                end

            else
                error('wrong number of arguments');
            end

        end

        %% Mixture Set methods

        function set.T(obj, temp)
            ct.impl.call('mMix_setTemperature', obj.mixID, temp);
        end

        function set.P(obj, pressure)
            ct.impl.call('mMix_setPressure', obj.mixID, pressure);
        end

        function setPhaseMoles(obj, n, moles)
            % Set the number of moles [kmol] of a phase in a mixture. ::
            %
            %     >> m.setPhaseMoles(n, moles)
            %
            % :param n:
            %     Phase number in the input.
            % :param moles:
            %     Number of moles to add.

            ct.impl.call('mMix_setPhaseMoles', obj.mixID, n - 1, moles);
        end

        function setSpeciesMoles(obj, moles)
            % Set the moles [kmol] of the species. ::
            %
            %     >> m.setSpeciesMoles(moles)
            %
            % Set the moles of multiple species. The moles may be specified either as a
            % string, or as an vector. If a vector is used, it must be dimensioned at
            % least as large as the total number of species in the mixture. Note that
            % the species may belong to any phase, and unspecified species are set to
            % zero. ::
            %
            %     >> mix.setSpeciesMoles('C(s):1.0, CH4:2.0, O2:0.2');
            %
            % :param moles:
            %     Vector or string specifying the moles of species.

            if isa(moles, 'double')
                l = length(moles);
                ct.impl.call('mMix_setMoles', obj.mixID, l, moles);
            elseif isa(moles, 'char')
                ct.impl.call('mMix_setMolesByName', obj.mixID, moles);
            else
                error('The input must be a vector or string!');
            end

        end

        function r = equilibrate(obj, XY, solver, rtol, maxsteps, maxiter, estimate_equil)
            % Set the mixture to a state of chemical equilibrium. ::
            %
            %     >> m.equilibrate(XY, solver, rtol, maxsteps, maxiter, estimate_equil)
            %
            % :param XY:
            %     Two-letter string specifying the two properties to hold
            %     fixed.  Currently, ``'TP'``, ``'HP'``, ``'TV'``, and ``'SP'`` are
            %     implemented. Default: ``'TP'``.
            % :param solver:
            %     Name of the solver to be used to equilibrate the phase.
            %     If solver = 'element_potential', the ChemEquil element potential
            %     solver will be used.
            %     If solver = 'vcs', the VCS solver will be used.
            %     If solver = 'gibbs', the MultiPhaseEquil solver will be used.
            %     If solver = 'auto', the solvers will be tried in order if the initial
            %     solver(s) fail.
            %     Default: ``'auto'``.
            % :param err:
            %     Error tolerance. Iteration will continue until :math:`\Delta\mu)/RT`
            %     is less than this value for each reaction. Default:
            %     1.0e-9. Note that this default is very conservative, and good
            %     equilibrium solutions may be obtained with larger error
            %     tolerances.
            % :param maxsteps:
            %     Maximum number of steps to take while solving the
            %     equilibrium problem for specified T and P. Default: 1000.
            % :param maxiter:
            %     Maximum number of temperature and/or pressure
            %     iterations.  This is only relevant if a property pair other
            %     than (T,P) is specified. Default: 200.
            % :param estimate_equil:
            %     For MultiPhaseEquil solver, an integer indicating whether the solver
            %     should estimate its own initial condition.
            %     If 0, the initial mole fraction vector in the ThermoPhase object
            %     is used as the initial condition.
            %     If 1, the initial mole fraction vector is used if the
            %     element abundances are satisfied.
            %     If -1, the initial mole fraction vector is thrown out,
            %     and an estimate is formulated.
            %     Default: 0.
            % :return:
            %     The error in the solution.

            arguments
                obj
                XY (1,1) string {mustBeMember(XY, ["TP", "HP", "TV", "SP"])} = "TP"
                solver (1,1) string {mustBeMember(solver, ["auto", ...
                                     "vcs", "gibbs", "element_potential"])} = "auto"
                rtol (1,1) double {mustBePositive} = 1.0e-9
                maxsteps (1,1) double {mustBeInteger, mustBePositive} = 1000
                maxiter (1,1) double {mustBeInteger, mustBePositive} = 200
                estimate_equil (1,1) double {mustBeInteger} = 0
            end

            r = ct.impl.call('mMix_equilibrate', obj.mixID, XY, solver, rtol, ...
                        maxsteps, maxiter, estimate_equil);
        end

    end

end
