classdef (Abstract) Domain1D < handle
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

        domainID = -1  % ID of the domain

        % The Solution object used to represent the contents of this domain.
        phase Solution

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

        function d = Domain1D(id)
            % Create a :mat:class:`Domain1D` object.
            arguments
                id (1,1) double {mustBeInteger}
            end
            d.domainID = id;
            phaseID = ctFunc('mDomain_phase', id);
            d.phase= Solution(phaseID);
        end

        %% Domain1D Class Destructor

        function delete(d)
            % Delete the :mat:class:`Domain1D` object.
            if d.domainID >= 0
                ctFunc('mDomain_del', d.domainID);
            end
        end

        %% Domain1D Utility Methods

        function updateState(d, loc)
            % Set state of associated phase to specified location.
            arguments
                d
                loc (1,1) double {mustBeInteger} = 0
            end
            ctFunc('mDomain_updateState', d.domainID, loc);
        end

        function info(d, rows, width)
            % Print a concise summary of a Domain.
            %
            %     >> d.info()
            %
            % :param rows:
            %       Maximum number of rendered rows; defaults to 10.
            % :param width:
            %       Maximum width of rendered output; default adjusts to terminal width.
            arguments
                d (1,1) Domain1D
                rows (1,1) double {mustBeInteger, mustBePositive} = 10
                width (1,1) double {mustBeInteger} = -1
            end
            if width < 0
                try
                    size = matlab.desktop.commandwindow.size;
                    width = size(1);
                catch ME
                    width = 100;
                end
            end
            disp(ctString('mDomain_info', d.domainID, rows, width));
            fprintf("\n")
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
            lower = ctFunc('mDomain_lowerBound', d.domainID, n);
            upper = ctFunc('mDomain_upperBound', d.domainID, n);
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
                n = ctFunc('mDomain_componentIndex', d.domainID, name);

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
                output = ctString('mDomain_componentName', d.domainID, id);
                s{i} = output;
            end

        end

        function i = get.domainIndex(d)
            i = ctFunc('mDomain_domainIndex', d.domainID) + 1;
        end

        function str = get.domainType(d)
            str = ctString('mDomain_type', d.domainID);
        end

        function n = get.nComponents(d)
            n = ctFunc('mDomain_nComponents', d.domainID);
        end

        function n = get.nPoints(d)
            n = ctFunc('mDomain_nPoints', d.domainID);
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
            rerr = ctFunc('mDomain_rtol', d.domainID, n);
            aerr = ctFunc('mDomain_atol', d.domainID, n);
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
            ctFunc('mDomain_setBounds', d.domainID, n - 1, lower, upper);
        end

        function setSteadyTolerances(d, rtol, atol, component)
            % Set the steady-state tolerances. ::
            %
            %     >>d.setSteadyTolerances(rtol, atol, component)
            %
            % :param d:
            %     Instance of class :mat:class:`Domain1D`.
            % :param rtol:
            %     Relative tolerance.
            % :param atol:
            %     Absolute tolerance.
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.

            if nargin < 4 | strcmp(component, 'default')
                ctFunc('mDomain_setSteadyTolerances', d.domainID, rtol, atol, -1);
            elseif iscell(component)
                for ii = 1:length(component)
                    n = d.componentIndex(component{ii});
                    ctFunc('mDomain_setSteadyTolerances', d.domainID, rtol, atol, n);
                end
            else
                n = d.componentIndex(component);
                ctFunc('mDomain_setSteadyTolerances', d.domainID, rtol, atol, n);
            end

        end

        function setTransientTolerances(d, rtol, atol, component)
            % Set the transient tolerances. ::
            %
            %     >> d.setTransientTolerances(rtol, atol, component)
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

            if nargin < 4 | strcmp(component, 'default')
                ctFunc('mDomain_setTransientTolerances', d.domainID, rtol, atol, -1);
            elseif iscell(component)
                for ii = 1:length(component)
                    n = d.componentIndex(component{ii});
                    ctFunc('mDomain_setTransientTolerances', d.domainID, rtol, atol, n);
                end
            else
                n = d.componentIndex(component);
                ctFunc('mDomain_setTransientTolerances', d.domainID, rtol, atol, n);
            end

        end

    end

end
