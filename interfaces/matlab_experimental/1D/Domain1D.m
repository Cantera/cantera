classdef Domain1D < handle

    properties
        dom_id
        domain_type
        T
    end

    methods
        %% Domain1D class constructor.

        function d = Domain1D(a, b, c)
            % parameter a:
            %    Integer type of domain. Possible values are:
            %    * 1 - Stagnation Flow
            %    * 2 - Inlet1D
            %    * 3 - Surf1D
            %    * 4 - Symm1D
            %    * 5 - Outlet1D
            %    * 6 - Reacting Surface
            %    * 8 - Sim1D
            %    * -2 - OutletRes
            %
            % parameter b:
            %    Instance of class 'Solution' (for a == 1) or 'Interface'
            %    (for a == 6). Not used for all other valid values of 'a'.
            %
            % parameter c:
            %    Integer, either 1 or 2, indicating whether an axisymmetric
            %    stagnation flow or a free flame should be created. If not
            %    specified, defaults to 1. Ignored if a!= 1.

            checklib;
            d.dom_id = -1;

            if nargin == 1
                if a == 2
                    d.dom_id = calllib(ct, 'inlet_new');
                elseif a == 3
                    d.dom_id = calllib(ct, 'surf_new');
                elseif a == 4
                    d.dom_id = calllib(ct, 'symm_new');
                elseif a == 5
                    d.dom_id = calllib(ct, 'outlet_new');
                elseif a == -2
                    d.dom_id = calllib(ct, 'outletres_new');
                else
                    error('Not enough arguments for that job number');
                end
            elseif nargin == 2
                % a stagnation flow
                if a == 1
                    if isa(b, 'Solution')
                        d.dom_id = calllib(ct, 'stflow_new', ...
                                           b.tp_id, b.kin_id, b.tr_id, 1);
                    else
                        error('Wrong argument type. Expecting instance of class Solution.');
                    end
                elseif a == 6
                    if isa(b, 'Interface')
                        d.dom_id = calllib(ct, 'reactingsurf_new');
                        calllib(ct, 'reactingsurf_setkineticsmgr', ...
                                d.dom_id, b.kin_id);
                    else
                        error('Wrong argument type. Expecting instance of class Interface.');
                    end
                else
                    error('Wrong object type.');
                end
            elseif nargin == 3
                if a == 1
                    if isa(b, 'Solution')
                        d.dom_id = calllib(ct, 'stflow_new', ...
                                           b.tp_id, b.kin_id, b.tr_id, c);
                    else
                        error('Wrong argument type. Expecting instance of class Solution.');
                    end
                else
                    error('Unknown domain type.');
                end
            end
%             if d.dom_id < 0
%                 error(geterr);
%             end
            d.domain_type = a;
        end

        %% Utility Methods

        function dom_clear(d)
            % Delete the Domain1D object
            checklib;
            calllib(ct, 'domain_del', d.dom_id);
        end

        %% Domain Methods

        function n = componentIndex(d, name)
            % Get the index of a component given its name
            %
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter name:
            %    String name of the component to look up. If a numeric
            %    value is passed, it will be returned directly.
            % return:
            %    Index of the component.

            checklib;
            if isa(name, 'double')
                n = name;
            else
                n = calllib(ct, 'domain_componentIndex', ...
                            d.dom_id, name);
                if n >= 0
                    n = n+1;
                end
            end
            if n <= 0
                error('Component not found');
            end
        end

        function n = componentName(d, index)
            % Get the name of a component given its index.
            %
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter n:
            %    Integer or vector of integers of components' names.
            % return:
            %    Cell array of component names.

            checklib;
            n = length(index);
            s = cell(m);
            for i = 1:n
                id = index(i)-1;
                buflen = calllib(ct, 'domain_componentName', ...
                                 d.dom_id, id, 0, 0);
                if buflen > 0
                    aa = char(zeros(1, buflen));
                    [out_buf, aa] = calllib(ct, 'domain_componentName', ...
                                            d.dom_id, id, buflen, aa);
                s{i} = aa;
            end
            end
        end

        function d = disableEnergy(d)
            % Disable the energy equation.
            %
            % parameter d:
            %    Instance of class 'Domain1D'.

            checklib;
            disp(' ');
            disp('Disabling the energy equation...');
            calllib(ct, 'stflow_solveEnergyEqn', d.dom_id, 0);
        end

        function i = domainIndex(d)
            % Get the domain index.
            % return:
            %    This function returns an integer flag denoting the
            %    location of the domain, beginning with 1 at the left.
            checklib;
            i = calllib(ct, 'domain_index', d.dom_id);
            if i >= 0
                i = i + 1;
            end
            if i <= 0
                error('Domain not found');
            end
        end

        function i = domainType(d)
            % Get the type of domain.
            % return:
            %    This function returns an integer flag denoting the domain
            %    type.
            checklib;
            i = calllib(ct, 'domain_type', d.dom_id);
        end

        function d = enableEnergy(d)
            % Enable the energy equation.
            % parameter d:
            %    Instance of class 'Domain1D'.

            checklib;
            disp(' ');
            disp('Enabling the energy equation...');
            calllib(ct, 'stflow_solveEnergyEqn', d.dom_id, 1);
        end

        function zz = gridPoints(d, n)
            % Get grid points from a domain.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter n:
            %    Optional, vector of grid points to be retrieved.
            % return:
            %    Vector of grid points.

            checklib;
            if nargin == 1
                np = d.nPoints;
                zz = zeros(1, np);
                for i = 1:np
                    zz(i) = calllib(ct, 'domain_grid', d.dom_id, i-1);
                end
            else
                m = length(n);
                zz = zeros(1, m);
                for i = 1:m
                    zz(i) = calllib(ct, 'domain_grid', d.dom_id, n(i)-1);
                end
            end
        end

        function a = isFlow(d)
            % Determine whether a domain is a flow.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % return:
            %    1 if the domain is a flow domain, and 0 otherwise.

            checklib;
            t = d.domainType;
            % See Domain1D.h for definitions of constants.
            if t < 100
                a = 1;
            else a = 0;
            end
        end

        function a = isInlet(d)
            % Determine whether a domain is an inlet.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % return:
            %    1 if the domain is an inlet, and 0 otherwise.

            checklib;
            t = d.domainType;
            % See Domain1D.h for definitions of constants.
            if t == 104
                a = 1;
            else a = 0;
            end
        end

        function a = isSurface(d)
            % Determine whether a domain is a surface.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % return:
            %    1 if the domain is a surface, and 0 otherwise.

            checklib;
            t = d.domainType;
            % See Domain1D.h for definitions of constants.
            if t == 102
                a = 1;
            else a = 0;
            end
        end

        function mdot = massFlux(d)
            % Get the mass flux.
            % return:
            %    The mass flux in the domain.
            checklib;
            mdot = calllib(ct, 'bdry_mdot', d.dom_id);
        end

        function y = massFraction(d, k)
            % Get the mass fraction of a species given its integer index.
            % This method returns the mass fraction of species 'k', where k
            % is the integer index of the species in the flow domain to
            % which the boundary domain is attached.
            %
            % parameter d:
            %    Instance of class 'Domain1D'
            % parameter k:
            %    Integer species index
            % return:
            %    Mass fraction of species

            checklib;
            if d.domainIndex == 0
                error('No flow domain attached!')
            end

            if d.isInlet
                y = calllib(ct, 'bdry_massFraction', d.dom_id, k-1);
            else error('Input domain must be an inlet');
            end
        end

        function n = nComponents(d)
            % Get the number of components.
            % return:
            %    Number of variables at each grid point.
            checklib;
            n = calllib(ct, 'domain_nComponents', d.dom_id);
        end

        function n = nPoints(d)
            % Get the number of grid points.
            % return:
            %    Integer number of grid points.
            checklib;
            n = calllib(ct, 'domain_nPoints', d.dom_id);
        end

        function setBounds(d, component, lower, upper)
            % Set bounds on the solution components.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter component:
            %    String, component to set the bounds on.
            % parameter lower:
            %    Lower bound.
            % parameter upper:
            %    Upper bound.
            checklib;
            n = d.componentIndex(component);
            calllib(ct, 'domain_setBounds', d.dom_id, n-1, lower, upper);
        end

        function setCoverageEqs(d, onoff)
            % Enable or disable solving the coverage equations.
            % parameter d:
            %    Instance of class 'Domain1D'
            % parameter onoff:
            %    string, one of 'on' or 'yes' to turn solving the coverage
            %    equations on. One of 'off' or 'no' to turn off the
            %    coverage equation.
            checklib;
            if d.domain_type ~= 6
                error('Wrong domain type. Expected a reacting surface domain.');
            end

            ion = -1;
            if isa(onoff,'char')
                if strcmp(onoff, 'on') || strcmp(onoff, 'yes')
                    ion = 1;
                elseif strcmp(onoff, 'off') || strcmp(onoff, 'no')
                    ion = 0;
                else
                    error(strcat('unknown option: ', onoff))
                end
            elseif isa(onoff, 'numeric')
                ion = onoff;
            end
            calllib(ct, 'reactingsurf_enableCoverageEqs', d.dom_id, ion);
        end

        function setFixedTempProfile(d, profile)
            % Set a fixed temperature profile to use when the energy
            % equation is not being solved. The profile must be entered as
            % an array of positions/temperatures, which may be in rows or
            % columns.
            %
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter profile:
            %    n x 2 or 2 x n array of 'n' points at which the
            %    temperature is specified.
            checklib;
            sz = size(profile);
            if sz(1) == 2
                l = length(profile(1, :));
                calllib(ct, 'stflow_setFixedTempProfile', d.dom_id, ...
                        l, profile(1, :), l, profile(2, :));
            elseif sz(2) == 2
                l = length(profile(:, 1));
                calllib(ct, 'stflow_setFixedTempProfile', d.dom_id, ...
                        l, profile(:, 1), l, profile(:, 2));
            else error('Wrong temperature profile array shape.');
            end
        end

        function setID(d, id)
            % Set the ID tag for a domain.
            % parameter id:
            %    String ID to assign.
            checklib;
            calllib(ct, 'domain_setID', d.dom_id, id);
        end

        function setMdot(d, mdot)
            % Set the mass flow rate.
            % parameter mdot:
            %    Mass flow rate.
            checklib;
            calllib(ct, 'bdry_setMdot', d.dom_id, mdot);
        end

        function setMoleFractions(d, x)
            % Set the mole fractions.
            % parameter x:
            %    String specifying the species and mole fractions in the
            %    format "'Spec:X,Spec2:X2'".
            checklib;
            calllib(ct, 'bdry_setMoleFractions', d.dom_id, x);
        end

        function setPressure(d, p)
            % Set the pressure.
            % parameter p:
            %    Pressure to be set. Unit: Pa.
            checklib;
            calllib(ct, 'stflow_setPressure', d.dom_id, p);
        end

        function setProfileD(d, n, p)
            % Set the profile of a component.
            % Convenience function to allow an instance of 'Domain1D' to
            % have a profile of its components set when it is part of a
            % 'Stack'.
            %
            % parameter d:
            %    Instance of class 'Domain1d'.
            % parameter n:
            %    Integer index of component, vector of component indices,
            %    string of component name, or cell array of strings of
            %    component names.
            % parameter p:
            %    n x 2 array, whose columns are the relative (normalized)
            %    positions and the component values at those points. The
            %    number of positions 'n' is arbitrary.

            checklib;
            if d.stack == 0
                error('Install domain in stack before calling setProfile.');
            end
            d.stack.setProfile(d.domainIndex, n, p);
        end

        function setSteadyTolerances(d, component, rtol, atol)
            % Set the steady-state tolerances.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter component:
            %    String or cell array of strings of component values whose
            %    tolerances should be set. If 'default' is specified, the
            %    tolerance of all components will be set.
            % parameter rtol:
            %    Relative tolerance.
            % parameter atol:
            %    Absolute tolerance.

            checklib;
            if strcmp(component, 'default')
                nc = d.nComponents;
                for ii = 1:nc
                    calllib(ct, 'domain_setSteadyTolerances', ...
                            d.dom_id, ii, rtol, atol);
                end
            elseif iscell(component)
                nc = length(component);
                for ii = 1:nc
                    n = d.componentIndex(component{ii});
                    calllib(ct, 'domain_setSteadyTolerances', ...
                            d.dom_id, n, rtol, atol);
                end
            else
                n = d.componentIndex(component);
                calllib(ct, 'domain_setSteadyTolerances', ...
                        d.dom_id, ii, rtol, atol);
            end
        end

        function setTransientTolerances(d, component, rtol, atol)
            % Set the transient tolerances.
            % parameter d:
            %    Instance of class 'Domain1D'.
            % parameter component:
            %    String or cell array of strings of component values whose
            %    tolerances should be set. If 'default' is specified, the
            %    tolerance of all components will be set.
            % parameter rtol:
            %    Relative tolerance.
            % parameter atol:
            %    Absolute tolerance.

            checklib;
            if strcmp(component, 'default')
                nc = d.nComponents;
                for ii = 1:nc
                    calllib(ct, 'domain_setTransientTolerances', ...
                            d.dom_id, ii, rtol, atol);
                end
            elseif iscell(component)
                nc = length(component);
                for ii = 1:nc
                    n = d.componentIndex(component{ii});
                    calllib(ct, 'domain_setTransientTolerances', ...
                            d.dom_id, n, rtol, atol);
                end
            else
                n = d.componentIndex(component);
                calllib(ct, 'domain_setTransientTolerances', ...
                        d.dom_id, ii, rtol, atol);
            end
        end

        function temperature = get.T(d)
            % Get the boundary temperature (K).
            checklib;
            temperature = calllib(ct, 'bdry_temperature', d.dom_id);
        end

        function set.T(d, t)
            % Set the temperature (K).
            checklib;
            if t <= 0
                error('The temperature must be positive');
            end
            calllib(ct, 'bdry_setTemperature', d.dom_id, t);
        end

        function setupGrid(d, grid)
            % Set up the solution grid.
            checklib;
            calllib(ct, 'domain_setupGrid', d.dom_id, numel(grid), grid);
        end

    end
end
