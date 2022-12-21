classdef Domain1D < handle
    % Domain1D Class
    %
    % d = Domain1D(a, b)
    %
    % :param a:
    %    String type of domain. Possible values are:
    %    `StagnationFlow`
    %    `AxisymmetricFlow`
    %    `Inlet1D`
    %    `Surf1D`
    %    `Symm1D`
    %    `Outlet1D`
    %    `ReactingSurface`
    %    `OutletRes`
    %
    % :param b:
    %     Instance of class :mat:class:`Solution` (for ``a == 1``)
    %     or :mat:class:`Interface` (for ``a == 6``). Not used for
    %     all other valid values of ``a``.
    %

    properties (SetAccess = immutable)

        domainID % ID of the domain

        type % Type of the domain
    end

    properties (SetAccess = public)

        T % Boundary Temperature. Units: K.

        P % Boundary Pressure. Units: Pa.
    end

    properties (SetAccess = protected)

        % Domain index.
        %
        % i = d.domainIndex
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     Integer flag denoting the location of the domain,
        %     beginning with 1 at the left.
        domainIndex
        %
        % Type of the domain.
        %
        % i = d.domainType
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     Integer flag denoting the domain type.
        domainType
        %
        % Determines whether a domain is a flow.
        %
        % a = d.isFlow
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     1 if the domain is a flow domain, and 0 otherwise.
        isFlow
        %
        % Determines whether a domain is an inlet.
        %
        % a = d.isInlet
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     1 if the domain is an inlet, and 0 otherwise.
        isInlet
        %
        % Determines if a domain is a surface.
        %
        % a = d.isSurface
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     1 if the domain is a surface, and 0 otherwise.
        isSurface
        %
        % Mass flux.
        %
        % mdot = d.massFlux
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     The mass flux in the domain.
        massFlux
        %
        % Number of components.
        %
        % n = d.nComponents
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     Number of variables at each grid point
        nComponents
        %
        % Get the number of grid points.
        %
        % n = d.nPoints
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`
        % :return:
        %     Integer number of grid points.
        nPoints

    end

    methods
        %% Domain1D Class Constructor.

        function d = Domain1D(a, b)
            checklib;

            if nargin == 1

                if strcmp(a, 'Inlet1D')
                    d.domainID = callct('inlet_new');
                elseif strcmp(a, 'Surf1D')
                    d.domainID = callct('surf_new');
                elseif strcmp(a, 'Symm1D')
                    d.domainID = callct('symm_new');
                elseif strcmp(a, 'Outlet1D')
                    d.domainID = callct('outlet_new');
                elseif strcmp(a, 'OutletRes')
                    d.domainID = callct('outletres_new');
                else
                    error('Not enough arguments for that job number');
                end

            elseif nargin == 2
                % a stagnation flow
                if strcmp(a, 'StagnationFlow')

                    if isa(b, 'Solution')
                        d.domainID = callct('stflow_new', ...
                                            b.tpID, b.kinID, b.trID, 1);
                    else
                        error('Wrong argument type. Expecting instance of class Solution.');
                    end

                elseif strcmp(a, 'AxisymmetricFlow')

                    if isa(b, 'Solution')
                        d.domainID = callct('stflow_new', ...
                                            b.tpID, b.kinID, b.trID, 2);
                    else
                        error('Wrong argument type. Expecting instance of class Solution.');
                    end

                elseif strcmp(a, 'ReactingSurface')

                    if isa(b, 'Interface')
                        d.domainID = callct('reactingsurf_new');
                        callct('reactingsurf_setkineticsmgr', ...
                                d.domainID, b.kinID);
                    else
                        error('Wrong argument type. Expecting instance of class Interface.');
                    end

                else
                    error('Wrong object type.');
                end

            end

            d.type = a;
        end

        %% Domain1D Class Destructor

        function delete(d)
            % Delete the C++ Domain1D object.

            callct('domain_del', d.domainID);
        end

        %% Domain1D Utility Methods

        function d = disableEnergy(d)
            % Disable the energy equation.
            %
            % d.disableEnergy
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            %
            disp(' ');
            disp('Disabling the energy equation...');
            callct('stflow_solveEnergyEqn', d.domainID, 0);
        end

        function d = enableEnergy(d)
            % Enable the energy equation.
            %
            % d.enableEnergy
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            %
            disp(' ');
            disp('Enabling the energy equation...');
            callct('stflow_solveEnergyEqn', d.domainID, 1);
        end

        function d = disableSoret(d)
            % Disable the diffusive mass fluxes due to the Soret effect.
            %
            % d.disableSoret
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            %
            disp(' ');
            disp('Disabling the Soret effect...');
            callct('stflow_enableSoret', d.domainID, 0);
        end

        function d = enableSoret(d)
            % Enable the diffusive mass fluxes due to the Soret effect.
            %
            % d.enableSoret
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            %
            disp(' ');
            disp('Disabling the Soret effect...');
            callct('stflow_enableSoret', d.domainID, 1);
        end

        %% Domain Get Methods

        function b = bounds(d, component)
            % Get the (lower, upper) bounds for a solution component.
            %
            % b = d.bounds(componoent)
            %
            % :param component:
            %    String name of the component for which the bounds are
            %    returned.
            % :return:
            %    1x2 Vector of the lower and upper bounds.

            n = d.componentIndex(component);
            lower = callct('domain_lowerBound', d.domainID, n);
            upper = callct('domain_upperBound', d.domainID, n);
            b = [lower, upper];
        end

        function n = componentIndex(d, name)
            %
            % Index of a component given its name.
            %
            % n = d.componentIndex(name)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param name:
            %     String name of the component to look up. If a numeric value
            %     is passed, it will be returned.
            % :return:
            %     Index of the component, or input numeric value.

            if isa(name, 'double')
                n = name;
            else
                n = callct('domain_componentIndex', ...
                    d.domainID, name);

                if n >= 0
                    n = n + 1;
                end

            end

            if n <= 0
                error('Component not found');
            end

        end

        function s = componentName(d, index)
            %
            % Name of a component given its index.
            %
            % n = d.componentName(index)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param index:
            %     Integer or vector of integers of component names
            %     to get.
            % :return:
            %     Cell array of component names.

            n = length(index);
            s = cell(1, n);

            for i = 1:n
                id = index(i) - 1;
                output = callct2('domain_componentName', d.domainID, id);
                s{i} = output;
            end

        end

        function i = get.domainIndex(d)
            i = callct('domain_index', d.domainID);

            if i >= 0
                i = i + 1;
            end

            if i <= 0
                error('Domain not found');
            end

        end

        function i = get.domainType(d)
            i = callct('domain_type', d.domainID);
        end

        function zz = gridPoints(d, n)
            %
            % Grid points from a domain.
            %
            % zz = d.gridPoints(n)
            %
            % :param d:
            %    Instance of class 'Domain1D'.
            % :param n:
            %    Optional, vector of grid points to be retrieved.
            % :return:
            %    Vector of grid points.

            if nargin == 1
                np = d.nPoints;
                zz = zeros(1, np);

                for i = 1:np
                    zz(i) = callct('domain_grid', d.domainID, i - 1);
                end

            else
                m = length(n);
                zz = zeros(1, m);

                for i = 1:m
                    zz(i) = callct('domain_grid', d.domainID, n(i) - 1);
                end

            end

        end

        function a = get.isFlow(d)
            t = d.domainType;
            % See Domain1D.h for definitions of constants.
            if t < 100
                a = 1;
            else a = 0;
            end

        end

        function a = get.isInlet(d)
            t = d.domainType;
            % See Domain1D.h for definitions of constants.
            if t == 104
                a = 1;
            else a = 0;
            end

        end

        function a = get.isSurface(d)
            t = d.domainType;
            % See Domain1D.h for definitions of constants.
            if t == 102
                a = 1;
            else a = 0;
            end

        end

        function mdot = get.massFlux(d)
            mdot = callct('bdry_mdot', d.domainID);
        end

        function y = massFraction(d, k)
            %
            % Get the mass fraction of a species given its integer index.
            %
            % y = d.massFraction(k)
            %
            % This method returns the mass fraction of species ``k``, where
            % k is the integer index of the species in the flow domain
            % to which the boundary domain is attached.
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param k:
            %     Integer species index
            % :return:
            %     Mass fraction of species

            if d.domainIndex == 0
                error('No flow domain attached!')
            end

            if d.isInlet
                y = callct('bdry_massFraction', d.domainID, k - 1);
            else error('Input domain must be an inlet');
            end

        end

        function n = get.nComponents(d)
            n = callct('domain_nComponents', d.domainID);
        end

        function n = get.nPoints(d)
            n = callct('domain_nPoints', d.domainID);
        end

        function tol = tolerances(d, component)
            %
            % Return the (relative, absolute) error tolerances for a
            % solution component.
            %
            % tol = d.tolerances(component)
            %
            % :param component:
            %    String name of the component for which the bounds are
            %    returned.
            % :return:
            %    1x2 Vector of the relative and absolute error tolerances.

            n = d.componentIndex(component);
            rerr = callct('domain_rtol', d.domainID, n);
            aerr = callct('domain_atol', d.domainID, n);
            tol = [rerr, aerr];
        end

        function temperature = get.T(d)
            temperature = callct('bdry_temperature', d.domainID);
        end

        function pressure = get.P(d)
            pressure = calllibt(ct, 'stflow_pressure', d.domainID);
        end

        %% Domain Set Methods

        function setBounds(d, component, lower, upper)
            % Set bounds on the solution components.
            %
            % d.setBounds(component, lower, upper)
            %
            % :param d:
            %    Instance of class 'Domain1D'.
            % :param component:
            %    String, component to set the bounds on.
            % :param lower:
            %    Lower bound.
            % :param upper:
            %    Upper bound.
            %

            n = d.componentIndex(component);
            callct('domain_setBounds', d.domainID, n - 1, lower, upper);
        end

        function setCoverageEqs(d, onoff)
            % Set bounds on the solution components.
            %
            % d.setBounds(component, lower, upper)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param component:
            %     String, component to set the bounds on
            % :param lower:
            %     Lower bound
            % :param upper:
            %     Upper bound
            %
            if ~strcmp(d.type, 'ReactingSurface')
                error('Wrong domain type. Expected a reacting surface domain.');
            end

            ion = -1;

            if isa(onoff, 'char')

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

            callct('reactingsurf_enableCoverageEqs', d.domainID, ion);
        end

        function setFixedTempProfile(d, profile)
            % Set a fixed temperature profile.
            %
            % d.setFixedTempProfile(profile)
            %
            % Set the temperature profile to use when the
            % energy equation is not being solved. The profile must be entered
            % as an array of positions / temperatures, which may be in rows or
            % columns.
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param profile:
            %     n x 2 or 2 x n array of ``n`` points at which the temperature
            %     is specified.
            %
            sz = size(profile);

            if sz(1) == 2
                l = length(profile(1, :));
                callct('stflow_setFixedTempProfile', d.domainID, ...
                        l, profile(1, :), l, profile(2, :));
            elseif sz(2) == 2
                l = length(profile(:, 1));
                callct('stflow_setFixedTempProfile', d.domainID, ...
                        l, profile(:, 1), l, profile(:, 2));
            else error('Wrong temperature profile array shape.');
            end

        end

        function setID(d, id)
            % Set the ID tag for a domain.
            %
            % d.setID(id)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param id:
            %     String ID to assign
            %
            callct('domain_setID', d.domainID, id);
        end

        function setMassFlowRate(d, mdot)
            % Set the mass flow rate.
            %
            % d.setMassFlowRate(mdot)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param mdot:
            %     Mass flow rate
            %
            callct('bdry_setMdot', d.domainID, mdot);
        end

        function setMoleFractions(d, x)
            % Set the mole fractions.
            %
            % d.setMoleFractions(x)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param x:
            %     String specifying the species and mole fractions in
            %     the format ``'SPEC:X,SPEC2:X2'``.
            %
            callct('bdry_setMoleFractions', d.domainID, x);
        end

        function setSteadyTolerances(d, component, rtol, atol)
            % Set the steady-state tolerances.
            %
            % d.setSteadyTolerances(component, rtol, atol)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.
            % :param rtol:
            %     Relative tolerance
            % :param atol:
            %     Absolute tolerance
            %
            if strcmp(component, 'default')
                nc = d.nComponents;

                for ii = 1:nc
                    callct('domain_setSteadyTolerances', ...
                            d.domainID, ii - 1, rtol, atol);
                end

            elseif iscell(component)
                nc = length(component);

                for ii = 1:nc
                    n = d.componentIndex(component{ii});
                    callct('domain_setSteadyTolerances', ...
                            d.domainID, n, rtol, atol);
                end

            else
                n = d.componentIndex(component);
                callct('domain_setSteadyTolerances', ...
                        d.domainID, n, rtol, atol);
            end

        end

        function setTransientTolerances(d, component, rtol, atol)
            % Set the transient tolerances.
            %
            % d.setTransientTolerances(component, rtol, atol)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.
            % :param rtol:
            %     Relative tolerance
            % :param atol:
            %     Absolute tolerance
            %
            if strcmp(component, 'default')
                nc = d.nComponents;

                for ii = 1:nc
                    callct('domain_setTransientTolerances', ...
                            d.domainID, ii - 1, rtol, atol);
                end

            elseif iscell(component)
                nc = length(component);

                for ii = 1:nc
                    n = d.componentIndex(component{ii});
                    callct('domain_setTransientTolerances', ...
                            d.domainID, n, rtol, atol);
                end

            else
                n = d.componentIndex(component);
                callct('domain_setTransientTolerances', ...
                        d.domainID, n, rtol, atol);
            end

        end

        function setTransport(d, itr)
            % Set the solution object used for calculating transport properties.
            %
            % d.setTransport(itr)
            %
            % :param itr:
            %    ID of the solution object for which transport properties
            %    are calculated.
            %
            callct('stflow_setTransport', d.domainID, itr);
        end

        function setupGrid(d, grid)
            % Set up the solution grid.
            %
            % d.setupGrid(grid)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`
            % :param grid:
            %
            callct('domain_setupGrid', d.domainID, numel(grid), grid);
        end

        function set.T(d, t)

            if t <= 0
                error('The temperature must be positive');
            end

            callct('bdry_setTemperature', d.domainID, t);
        end

        function set.P(d, p)

            if p <= 0
                error('The pressure must be positive');
            end

            callct('stflow_setPressure', d.domainID, p);
        end

    end

end
