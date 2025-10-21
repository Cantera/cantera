classdef Mixture < handle
    % Mixture Class ::
    %
    %     >> m = Mixture(phases)
    %
    % Class :mat:class:`Mixture` represents mixtures of one or more phases of matter.
    % To construct a mixture, supply a cell array of phases and mole numbers ::
    %
    %     >> gas = Solution('gas.yaml');
    %     >> graphite = Solution('graphite.yaml');
    %     >> mix = Mixture({gas, 1.0; graphite, 0.1});
    %
    % Phases may also be added later using the addPhase method ::
    %
    %     >> water = Solution('water.yaml');
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
    % :return:
    %     Instance of class :mat:class:`Mixture`.
    %

    properties (SetAccess = immutable)
        mixID
    end

    properties (SetAccess = public)
        T % Temperature. Units: K.
        P % Pressure. Units: Pa.
    end

    properties (SetAccess = protected)

        phases % Phases in the mixture

        nElements % Number of elements in a mixture.

        nPhases % Number of phases in a mixture.

        nSpecies % Number of species in a mixture.

        chemPotentials % Chemical potentials of species in the mixture. Units: J/kmol.
    end

    methods
        %% Mixture Class Constructor

        function m = Mixture(phases)
            % Create a :mat:class:`Mixture` object.

            ctIsLoaded;

            if nargin > 1
                error('Mixture: wrong number of arguments');
            end

            % Create an empty mixture.
            m.mixID = ctFunc('mix_new');
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

        function delete(m)
            % Delete the :mat:class:`Mixture` object.
            ctFunc('mix_del', m.mixID);
        end

        %% Mixture Utility methods

        function display(m)
            % Display the state of the mixture on the terminal.

            ctFunc('mix_updatePhases', m.mixID);
            [np, nc] = size(m.phases);

            for n = 1:np
                s = [sprintf('\n*******************    Phase %d', n) ...
                     sprintf('    ******************************\n\n Moles: %12.6g', ...
                     m.phaseMoles(n))];
                disp(s);
                display(m.phases{n, 1});
            end

        end

        function addPhase(m, phase, moles)
            % Add a phase to a mixture. ::
            %
            %     >> m.addPhase(phase, moles)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture` to which phases
            %     should be added.
            % :param phase:
            %     Instance of class :mat:class:`ThermoPhase` which should be added.
            % :param moles:
            %     Number of moles of the ``phase`` to be added to this mixture.
            %     Units: kmol.

            if ~isa(phase, 'ThermoPhase')
                error('Phase object of wrong type.');
            end

            if ~isa(moles, 'numeric')
                error('Number of moles must be numeric');
            end

            if moles < 0.0
                error('Negative moles');
            end

            ctFunc('mix_addPhase', m.mixID, phase.tpID, moles);

        end

        %% Mixture Get methods

        function temperature = get.T(m)
            temperature = ctFunc('mix_temperature', m.mixID);
        end

        function pressure = get.P(m)
            pressure = ctFunc('mix_pressure', m.mixID);
        end

        function n = get.nElements(m)
            n = ctFunc('mix_nElements', m.mixID);
        end

        function n = get.nPhases(m)
            n = ctFunc('mix_nPhases', m.mixID);
        end

        function n = get.nSpecies(m)
            n = ctFunc('mix_nSpecies', m.mixID);
        end

        function mu = get.chemPotentials(m)
            nsp = m.nSpecies;
            mu = ctArray('mix_getChemPotentials', nsp, m.mixID);
        end

        function n = nAtoms(m, e)
            % Number of atoms of an element in a mixture. ::
            %
            %     >> n = m.nAtoms(e)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param e:
            %     Index of element.
            % :return:
            %     Number of atoms for element e.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.

            n = ctFunc('mix_nPhases', m.mixID, k - 1, e - 1);
        end

        function n = elementIndex(m, name)
            % Index of an element. ::
            %
            %     >> n = m.elementIndex(name)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param name:
            %     Name of the element whose index is desired.
            % :return:
            %     Index of element with name ``name``.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.

            n = ctFunc('mix_elementIndex', m.mixID, name) + 1;
        end

        function n = speciesIndex(m, k, p)
            % Index of a species in a mixture. ::
            %
            %     >> n = m.speciesIndex(k, p)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param name:
            %     Name of the species whose index is desired.
            % :return:
            %     Index of species with name ``name``.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.

            n = ctFunc('mix_speciesIndex', m.mixID, k - 1, p - 1) + 1;
        end

        function moles = elementMoles(m, e)
            % Number of moles of an element in a mixture. ::
            %
            %     >> moles = m.elementMoles(e)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param e:
            %    Integer element number.
            % :return:
            %    Moles of element number 'e'. If input 'e' is empty, return
            %    moles of every element in the mixture. Unit: kmol.

            if nargin == 2
                moles = ctFunc('mix_elementMoles', m.mixID, e);
            elseif nargin == 1
                nel = m.nElements;
                moles = zeros(1, nel);

                for i = 1:nel
                    moles(i) = ctFunc('mix_elementMoles', m.mixID, i-1);
                end

            else
                error('wrong number of arguments');
            end

        end

        function moles = phaseMoles(m, n)
            % Number of moles of a phase in a mixture. ::
            %
            %     >> moles = m.phaseMoles(n)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param n:
            %    Integer phase number.
            % :return:
            %    Moles of element number 'n'. If input 'n' is empty, return
            %    moles of every element in the mixture. Unit: kmol.

            if nargin == 2
                moles = ctFunc('mix_phaseMoles', m.mixID, n);
            elseif nargin == 1
                np = m.nPhases;
                moles = zeros(1, np);

                for i = 1:np
                    moles(i) = ctFunc('mix_phaseMoles', m.mixID, i-1);
                end

            else
                error('wrong number of arguments');
            end

        end

        function moles = speciesMoles(m, k)
            % Number of moles of a species in a mixture. ::
            %
            %     >> moles = m.speciesMoles(k)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param k:
            %    Integer species number.
            % :return:
            %    Moles of species number 'k'. If input 'k' is empty, return
            %    moles of every species in the mixture. Unit: kmol.

            if nargin == 2
                moles = ctFunc('mix_speciesMoles', m.mixID, k);
            elseif nargin == 1
                nsp = m.nSpecies;
                moles = zeros(1, nsp);

                for i = 1:nsp
                    moles(i) = ctFunc('mix_speciesMoles', m.mixID, i-1);
                end

            else
                error('wrong number of arguments');
            end

        end

        %% Mixture Set methods

        function set.T(m, temp)
            ctFunc('mix_setTemperature', m.mixID, temp);
        end

        function set.P(m, pressure)
            ctFunc('mix_setPressure', m.mixID, pressure);
        end

        function setPhaseMoles(m, n, moles)
            % Set the number of moles of a phase in a mixture. ::
            %
            %     >> m.setPhaseMoles(n, moles)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param n:
            %     Phase number in the input.
            % :param moles:
            %     Number of moles to add. Units: kmol.

            ctFunc('mix_setPhaseMoles', m.mixID, n - 1, moles);
        end

        function setSpeciesMoles(m, moles)
            % Set the moles of the species. ::
            %
            %     >> m.setSpeciesMoles(moles)
            %
            % Set the moles of the species in kmol. The moles may be specified
            % either as a string, or as an vector. If a vector is used,
            % it must be dimensioned at least as large as the total number
            % of species in the mixture. Note that the species may belong to any
            % phase, and unspecified species are set to zero. ::
            %
            %     >> mix.setSpeciesMoles('C(s):1.0, CH4:2.0, O2:0.2');
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
            % :param moles:
            %     Vector or string specifying the moles of species.

            if isa(moles, 'double')
                l = length(moles);
                ctFunc('mix_setMoles', m.mixID, l, moles);
            elseif isa(moles, 'char')
                ctFunc('mix_setMolesByName', m.mixID, moles);
            else
                error('The input must be a vector or string!');
            end

        end

        function r = equilibrate(m, XY, solver, rtol, maxsteps, maxiter, estimate_equil)
            % Set the mixture to a state of chemical equilibrium. ::
            %
            %     >> m.equilibrate(XY, solver, rtol, maxsteps, maxiter, estimate_equil)
            %
            % :param m:
            %     Instance of class :mat:class:`Mixture`.
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
                m
                XY (1,1) string {mustBeMember(XY, ["TP", "HP", "TV", "SP"])} = "TP"
                solver (1,1) string {mustBeMember(solver, ["auto", ...
                                     "vcs", "gibbs", "element_potential"])} = "auto"
                rtol (1,1) double {mustBePositive} = 1.0e-9
                maxsteps (1,1) double {mustBeInteger, mustBePositive} = 1000
                maxiter (1,1) double {mustBeInteger, mustBePositive} = 200
                estimate_equil (1,1) double {mustBeInteger} = 0
            end

            r = ctFunc('mix_equilibrate', m.mixID, XY, solver, rtol, ...
                        maxsteps, maxiter, estimate_equil);
        end

    end

end
