classdef Mixture < handle

    properties
        mixindex
        phases
        T
        P
    end

    methods
        %% Mixture class constructor

        function m = Mixture(phases)
            % To construct a mixture, supply a cell array of phases and mole
            % numbers:
            %    >> gas = Solution('gas.cti');
            %    >> graphite = Solution('graphite.cti');
            %    >> mix = Mixture({gas, 1.0; graphite, 0.1});
            %
            % Phases may also be added later using the addPhase method:
            %    >> water = Solution('water.cti');
            %    >> addPhase(mix, water, 3.0);
            %
            % Note that the objects representing each phase compute only
            % the intensive state of the phase - they do not store any
            % information on the amount of this phase. Mixture objects, on
            % the other hand, represent full extensive state.
            %
            % Mixture objects are 'lightweight'in the sense that they do
            % not store parameters needed to compute thermodynamic or
            % kinetic properties of the phases. These are contained in the
            % ('heavyweight') phase objects. Multiple mixture objects are
            % constructed using the same set of phase objects. Each one
            % store its own state information locally, and syncrhonizes the
            % phase objects whenever itrequires phase properties.
            %
            % parameter phases:
            %    Cell array of phases and mole numbers.
            % return:
            %    Instance of class 'Mixture'.

            checklib;

            if nargin > 1
                error('Mixture: wrong number of arguments');
            end

            % Create an empty mixture.
            m.mixindex = calllib(ct, 'mix_new');
            m.phases = phases;

            % If phases are supplied, add them
            if nargin == 1
                if ~isa(phases, 'cell')
                    error('Enter phasesas a cell array.');
                end

                % First column contains the phase objects, and the second
                % column contains the mole numbers of each phase.
                [np, nc] = size(phases);
                if nc ~= 2
                    error('Cell array of phases should have each phase on a new row');
                end
                for n = 1:np
                    m.addPhase(phases{n, 1}, phases{n, 2});
                end
                m.T = phases{n, 1}.T;
                m.P = phases{n, 1}.P;
            end
        end

        %% Utility methods

        function display(m)
            % Display the state of the mixture on the terminal.

            calllib(ct, 'mix_updatePhases', m.mixindex);
            [np, nc] = size(m.phases);
            for n = 1:np
                s = [sprintf('\n*******************    Phase %d', n) ...
                    sprintf('    ******************************\n\n Moles: %12.6g', ...
                            phaseMoles(m, n))];
                disp(s);
                display(m.phases{n, 1});
            end
        end

        function clear(m)
            % Delete the MultiPhase object.

            checklib;
            calllib(ct, 'mix_del', m.mixindex);
        end

        %% Mixture Get methods

        function addPhase(m, phase, moles)
            % Add a phase to the mixture
            %
            % parameter m:
            %    Instance of class 'Mixture' to which phases is added.
            % parameter phase:
            %    Instance of class 'ThermoPhase' which should be added.
            % parameter moles:
            %    Number of moles of the phase to be added. Unit: kmol.

            checklib;

            if ~isa(phase, 'ThermoPhase')
                error('Phase object of wrong type.');
            end
            if ~isa(moles, 'numeric')
                error('Number of moles must be numeric');
            end
            if moles < 0.0
                error('Negative moles');
            end

            iok = calllib(ct, 'mix_addPhase', m.mixindex, phase.tp_id, ...
                          moles);
            if iok < 0
                error('Error adding phase');
            end
        end

        function temperature = get.T(m)
            % Get the temperature of the mixture.
            %
            % return:
            %    Temperature in K.

            checklib;
            temperature = calllib(ct, 'mix_temperature', m.mixindex);
        end

        function pressure = get.P(m)
            % Get the pressure of themixture.
            %
            % return:
            %    Pressure in Pa.

            checklib;
            pressure = calllib(ct, 'mix_pressure', m.mixindex);
        end

        function n = nPhases(m)
            % Get the number of phases in the mixture.
            checklib;
            n = calllib(ct, 'mix_nPhases', m.mixindex);
        end

        function n = nElements(m)
            % Get the number of elements in the mixture.
            checklib;
            n = calllib(ct, 'mix_nElements', m.mixindex);
        end

        function n = nSpecies(m)
            % Get the number of species in the mixture.
            checklib;
            n = calllib(ct, 'mix_nSpecies', m.mixindex);
        end

        function n = elementIndex(m, name)
            % Get the index of element 'name'.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.
            checklib;
            n = calllib(ct, 'mix_elementIndex', m.mixindex, name) + 1;
        end

        function n = speciesIndex(m, k, p)
            % Get the index of element 'name'.
            %
            % Note: In keeping with the conventions used by Matlab, the
            % indices start from 1 instead of 0 as in Cantera C++ and
            % Python interfaces.
            checklib;
            n = calllib(ct, 'mix_speciesIndex', m.mixindex, k-1, p-1) + 1;
            % check back on this one!
        end

        function moles = phaseMoles(m, n)
            % Get the number of moles of a phase in a mixture.
            %
            % parameter n:
            %    Integer phase number in the input.
            % return:
            %    Moles of phase number 'n'. Unit: kmol.

            checklib;
            if nargin == 2
                moles = calllib(ct, 'mix_phaseMoles', m.mixindex, n);
            elseif nargin == 1
                np = m.nPhases;
                moles = zeros(1, np);
                for i = 1:np
                    moles(i) = calllib(ct, 'mix_phaseMoles', ...
                                       m.mixindex, i);
                end
            else error('wrong number of arguments');
            end
        end

        function mu = chemPotentials(m)
            % Get the chemical potentials of species in the mixture.
            %
            % return:
            %    Vector of chemical potentials. Unit: J/kmol.

            checklib;
            nsp = m.nSpecies;
            xx = zeros(1, nsp);
            ptr = libpointer('doublePtr', xx);
            calllib(ct, 'mix_getChemPotential', m.mixindex, nsp, ptr);
            mu = ptr.Value;
        end

        %% Mixture Set methods

        function m = set.T(m, temp)
            % Set the mixture temperature.
            %
            % parameter temp:
            %    Temperature to set. Unit: K.
            checklib;
            calllib(ct, 'mix_setTemperature', m.mixindex, temp);
        end

        function m = set.P(m, pressure)
            % Set the mixture pressure.
            %
            % parameter pressure:
            %    Pressure to set. Unit: Pa.
            checklib;
            calllib(ct, 'mix_setPressure', m.mixindex, pressure);
        end

        function setPhaseMoles(m, n, moles)
            % Set the number of moles of phase n in the mixture.
            %
            % parameter n:
            %    Phase number.
            % parameter moles:
            %    Number of moles to set. Unit: kmol.
            checklib;
            calllib(ct, 'mix_setPhaseMoles', m.mixindex, n-1, moles);
        end

        function setSpeciesMoles(m, moles)
            % Set the moles of the species. The moles may be specified as a
            % string or as a vector. If a vector is used, it must be
            % dimensioned at least as large as the total number of species
            % in the mixture. Note that the species may belong to any
            % phase, and unspecified species are set to zero.
            %
            % parameter moles:
            %    Vector or string specifying the moles of species.
            checklib;
            calllib(ct, 'mix_setMolesByName', m.mixindex, moles);
            % check back on this one!
        end

        function r = equilibrate(m, XY, err, maxsteps, maxiter, loglevel)
            % Set the mixture to a state of chemical equilibrium.
            %
            % This method uses a version of the VCS algorithm to find the
            % composition that minimizes the total Gibbs free energy of the
            % mixture, subject to element conservation constraints. For a
            % description of the theory, see Smith and Missen, "Chemical
            % Reaction Equilibrium". The VCS algorithm is implemented in
            % Cantera kernel class MultiPhaseEquil.
            %
            % The VCS algorithm solves for the equilibrium composition for
            % specified temperature and pressure. If any other property
            % pair other than 'TP' is specified, then an outer iteration
            % loop is used to adjust T and/or P so that the specified
            % property values are obtained:
            %     >> equilibrate(mix, 'TP);
            %     >> equilibrate('TP', 1.0e-6, 500);
            %
            % parameter XY:
            %    Two-letter string specifying the two properties to hold
            %    fixed. Currently 'TP', 'HP', 'TV', and 'SP' have been
            %    implemented. Default: 'TP'.
            % parameter err:
            %    Error tolerance. Iteration will continue until delta_Mu/RT
            %    is less than this value for each reaction. Default:
            %    1.0e-9.
            % parameter maxsteps:
            %    Maximum number of steps to take while solving the
            %    equilibrium problem for specified T and P. Default: 1000.
            % parameter maxiter:
            %    Maximum number of temperature and/or pressure iterations.
            %    This is only relevant if a property pair other than (T,
            %    P)is specified. Default: 200.
            % parameter loglevel:
            %    Set to a value > 0 to write diagnostic output. Larger
            %    values generate more detailed information.
            % return:
            %    The error in the solution.

            checklib;

            if nargin < 6
                loglevel = 0;
            end
            if nargin < 5
                maxiter = 200;
            end
            if nargin < 4
                maxsteps = 1000;
            end
            if nargin < 3
                err = 1.0e-9;
            end
            if nargin < 2
                XY = 'TP'
            end
            r = calllib(ct, 'mix_equilibrate', m.mixindex, XY, err, ...
                        maxsteps, maxiter, loglevel);
        end

    end
end
