classdef Domain1D < handle
    % Domain1D Class ::
    %
    %     >> d = Domain1D(type, phase, id)
    %
    % :param type:
    %    String type of domain. Possible values are:
    %      - `axisymmetric-flow`
    %      - `free-flow`
    %      - `inlet`
    %      - `outlet`
    %      - `reacting-surface`
    %      - `surface`
    %      - `symmetry-plane`
    %      - `outlet-reservoir`
    % :param phase:
    %     Instance of :mat:class:`Solution` or :mat:class:`Interface`.
    % :param id:
    %     String ID of the domain.
    % :return:
    %     Instance of class :mat:class:`Domain1D`

    properties (SetAccess = immutable)

        domainID % ID of the domain

    end

    properties (SetAccess = public)

        % Boolean flag indicating whether the energy equation is enabled.
        energyEnabled

        % Boolean flag indicating whether the diffusive mass fluxes due to the Soret
        % effect is enabled.
        soretEnabled

        % ID of the solution object used for calculating transport properties.
        transport
    end

    properties (SetAccess = protected)

        % Domain index. ::
        %
        %     >> i = d.domainIndex
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`.
        % :return:
        %     Integer flag denoting the location of the domain,
        %     beginning with 1 at the left.
        domainIndex

        % Type of the domain. ::
        %
        %     >> i = d.domainType
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`.
        % :return:
        %     Integer flag denoting the domain type.
        domainType

        % Number of components. ::
        %
        %     >> n = d.nComponents
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`.
        % :return:
        %     Number of variables at each grid point.
        nComponents

        % Get the number of grid points. ::
        %
        %     >> n = d.nPoints
        %
        % :param d:
        %     Instance of class :mat:class:`Domain1D`.
        % :return:
        %     Integer number of grid points.
        nPoints

    end

    methods
        %% Domain1D Class Constructor.

        function d = Domain1D(type, phase, id)
            % Create a :mat:class:`Domain1D` object.

            ctIsLoaded;

            d.domainID = ctFunc('domain_new', type, phase.solnID, id);

        end

        %% Domain1D Class Destructor

        function delete(d)
            % Delete the :mat:class:`Domain1D` object.

            if isempty(d.domainID)
                return
            end
            ctFunc('domain_del', d.domainID);
        end

        %% Domain1D Utility Methods

        function set.energyEnabled(d, flag)
            d.energyEnabled = flag;
            ctFunc('flow1D_solveEnergyEqn', d.domainID, int8(flag));
        end

        function set.soretEnabled(d, flag)
            d.soretEnabled = flag;
            ctFunc('flow1D_enableSoret', d.domainID, int8(flag));
        end

        %% Domain Get Methods

        function b = bounds(d, component)
            % Get the (lower, upper) bounds for a solution component. ::
            %
            %     >> b = d.bounds(component)
            %
            % :param component:
            %    String name of the component for which the bounds are returned.
            % :return:
            %    :math:`1\times 2` vector of the lower and upper bounds.

            n = d.componentIndex(component);
            lower = ctFunc('domain_lowerBound', d.domainID, n);
            upper = ctFunc('domain_upperBound', d.domainID, n);
            b = [lower, upper];
        end

        function n = componentIndex(d, name)
            % Index of a component given its name. ::
            %
            %     >>n = d.componentIndex(name)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param name:
            %     String name of the component to look up. If a numeric value
            %     is passed, it will be returned.
            % :return:
            %     Index of the component, or input numeric value.

            if isa(name, 'double')
                n = name;
            else
                n = ctFunc('domain_componentIndex', d.domainID, name);

                if n >= 0
                    n = n + 1;
                end

            end

            if n <= 0
                error('Component not found');
            end

        end

        function s = componentName(d, index)
            % Name of a component given its index. ::
            %
            %     >> n = d.componentName(index)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param index:
            %     Integer or vector of integers of component names to get.
            % :return:
            %     Cell array of component names.

            n = length(index);
            s = cell(1, n);

            for i = 1:n
                id = index(i) - 1;
                output = ctString('domain_componentName', d.domainID, id);
                s{i} = output;
            end

        end

        function i = get.domainIndex(d)
            i = ctFunc('domain_index', d.domainID) + 1;

            if i <= 0
                error('Domain not found');
            end

        end

        function str = get.domainType(d)
            str = ctString('domain_type', d.domainID);
        end

        function zz = gridPoints(d, n)
            % Grid points from a domain. ::
            %
            %     >> zz = d.gridPoints(n)
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
                    zz(i) = ctFunc('domain_grid', d.domainID, i - 1);
                end

            else
                m = length(n);
                zz = zeros(1, m);

                for i = 1:m
                    zz(i) = ctFunc('domain_grid', d.domainID, n(i) - 1);
                end

            end

        end

        function n = get.nComponents(d)
            n = ctFunc('domain_nComponents', d.domainID);
        end

        function n = get.nPoints(d)
            n = ctFunc('domain_nPoints', d.domainID);
        end

        function tol = tolerances(d, component)
            % Return the (relative, absolute) error tolerances for a
            % solution component. ::
            %
            %     >> tol = d.tolerances(component)
            %
            % :param component:
            %    String name of the component for which the bounds are returned.
            % :return:
            %    :math:`1\times 2` vector of the relative and absolute error tolerances.

            n = d.componentIndex(component);
            rerr = ctFunc('domain_rtol', d.domainID, n);
            aerr = ctFunc('domain_atol', d.domainID, n);
            tol = [rerr, aerr];
        end

        %% Domain Set Methods

        function setBounds(d, component, lower, upper)
            % Set bounds on the solution components. ::
            %
            %     >> d.setBounds(component, lower, upper)
            %
            % :param d:
            %    Instance of class :mat:class:`Domain1D`.
            % :param component:
            %    String, component to set the bounds on.
            % :param lower:
            %    Lower bound.
            % :param upper:

            n = d.componentIndex(component);
            ctFunc('domain_setBounds', d.domainID, n - 1, lower, upper);
        end



        function setSteadyTolerances(d, component, rtol, atol)
            % Set the steady-state tolerances. ::
            %
            %     >>d.setSteadyTolerances(component, rtol, atol)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.
            % :param rtol:
            %     Relative tolerance.
            % :param atol:
            %     Absolute tolerance.

            if strcmp(component, 'default')
                nc = d.nComponents;

                for ii = 1:nc
                    ctFunc('domain_setSteadyTolerances', ...
                            d.domainID, ii - 1, rtol, atol);
                end

            elseif iscell(component)
                nc = length(component);

                for ii = 1:nc
                    n = d.componentIndex(component{ii});
                    ctFunc('domain_setSteadyTolerances', d.domainID, n, rtol, atol);
                end

            else
                n = d.componentIndex(component);
                ctFunc('domain_setSteadyTolerances', d.domainID, n, rtol, atol);
            end

        end

        function setTransientTolerances(d, component, rtol, atol)
            % Set the transient tolerances. ::
            %
            %     >> d.setTransientTolerances(component, rtol, atol)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.
            % :param rtol:
            %     Relative tolerance.
            % :param atol:
            %     Absolute tolerance.

            if strcmp(component, 'default')
                nc = d.nComponents;

                for ii = 1:nc
                    ctFunc('domain_setTransientTolerances', ...
                            d.domainID, ii - 1, rtol, atol);
                end

            elseif iscell(component)
                nc = length(component);

                for ii = 1:nc
                    n = d.componentIndex(component{ii});
                    ctFunc('domain_setTransientTolerances', ...
                            d.domainID, n, rtol, atol);
                end

            else
                n = d.componentIndex(component);
                ctFunc('domain_setTransientTolerances', ...
                        d.domainID, n, rtol, atol);
            end

        end

        function set.transport(d, itr)
            ctFunc('flow1D_setTransport', d.domainID, itr);
        end

        function setupGrid(d, grid)
            % Set up the solution grid. ::
            %
            %     >> d.setupGrid(grid)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param grid:
            %     The grid for this domain.
            ctFunc('domain_setupGrid', d.domainID, numel(grid), grid);
        end

    end

end
