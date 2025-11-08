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

    properties (SetAccess = immutable)

        domainID = -1  % ID of the domain

        % The Solution object used to represent the contents of this domain.
        phase ct.Solution

    end

    properties (SetAccess = protected)

        % Integer denoting the location of the domain, beginning with 1 at the left.
        domainIndex

        % String denoting the domain type.
        domainType

        % Number of state components at each grid point.
        nComponents

        % Number of grid points in the domain.
        nPoints

    end

    methods
        %% Domain1D Class Constructor.

        function d = Domain1D(id)
            arguments
                id (1,1) double {mustBeInteger}
            end
            d.domainID = id;
            phaseID = ct.impl.call('mDomain_phase', id);
            d.phase = ct.Solution(phaseID);
        end

        %% Domain1D Class Destructor

        function delete(obj)
            % Delete the :mat:class:`Domain1D` object.
            if obj.domainID >= 0
                ct.impl.call('mDomain_del', obj.domainID);
            end
        end

        %% Domain1D Utility Methods

        function updateState(obj, loc)
            % Set state of associated phase to specified location.
            arguments
                obj
                loc (1,1) double {mustBeInteger} = 0
            end
            ct.impl.call('mDomain_updateState', obj.domainID, loc);
        end

        function info(obj, rows, width)
            % Print a concise summary of a Domain.
            %
            %     >> d.info()
            %
            % :param rows:
            %       Maximum number of rendered rows; defaults to 10.
            % :param width:
            %       Maximum width of rendered output; default adjusts to terminal width.
            arguments
                obj (1,1) ct.Domain1D
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
            disp(ct.impl.getString('mDomain_info', obj.domainID, rows, width));
            fprintf("\n")
        end

        %% Domain Get Methods

        function b = bounds(obj, component)
            % Get the (lower, upper) bounds for a solution component. ::
            %
            %     >> b = d.bounds(component)
            %
            % :param component:
            %    String name of the component for which the bounds are returned.
            % :return:
            %    :math:`1\times 2` vector of the lower and upper bounds.

            n = obj.componentIndex(component);
            lower = ct.impl.call('mDomain_lowerBound', obj.domainID, n);
            upper = ct.impl.call('mDomain_upperBound', obj.domainID, n);
            b = [lower, upper];
        end

        function n = componentIndex(obj, name)
            % Index of a component given its name. ::
            %
            %     >>n = d.componentIndex(name)
            %
            % :param name:
            %     String name of the component to look up. If a numeric value
            %     is passed, it will be returned.
            % :return:
            %     Index of the component, or input numeric value.

            if isa(name, 'double')
                n = name;
            else
                n = ct.impl.call('mDomain_componentIndex', obj.domainID, name);

                if n >= 0
                    n = n + 1;
                end

            end

            if n <= 0
                error('Component not found');
            end

        end

        function s = componentName(obj, index)
            % Name of a component given its index. ::
            %
            %     >> n = d.componentName(index)
            %
            % :param index:
            %     Integer or vector of integers of component names to get.
            % :return:
            %     Cell array of component names.

            n = length(index);
            s = cell(1, n);

            for i = 1:n
                id = index(i) - 1;
                output = ct.impl.getString('mDomain_componentName', obj.domainID, id);
                s{i} = output;
            end

        end

        function i = get.domainIndex(obj)
            i = ct.impl.call('mDomain_domainIndex', obj.domainID) + 1;
        end

        function str = get.domainType(obj)
            str = ct.impl.getString('mDomain_type', obj.domainID);
        end

        function n = get.nComponents(obj)
            n = ct.impl.call('mDomain_nComponents', obj.domainID);
        end

        function n = get.nPoints(obj)
            n = ct.impl.call('mDomain_nPoints', obj.domainID);
        end

        function tol = tolerances(obj, component)
            % Return the (relative, absolute) error tolerances for a
            % solution component. ::
            %
            %     >> tol = d.tolerances(component)
            %
            % :param component:
            %    String name of the component for which the bounds are returned.
            % :return:
            %    :math:`1\times 2` vector of the relative and absolute error tolerances.

            n = obj.componentIndex(component);
            rerr = ct.impl.call('mDomain_rtol', obj.domainID, n);
            aerr = ct.impl.call('mDomain_atol', obj.domainID, n);
            tol = [rerr, aerr];
        end

        %% Domain Set Methods

        function setBounds(obj, component, lower, upper)
            % Set bounds on the solution components. ::
            %
            %     >> d.setBounds(component, lower, upper)
            %
            % :param component:
            %    String, component to set the bounds on.
            % :param lower:
            %    Lower bound.
            % :param upper:

            n = obj.componentIndex(component);
            ct.impl.call('mDomain_setBounds', obj.domainID, n - 1, lower, upper);
        end

        function setSteadyTolerances(obj, rtol, atol, component)
            % Set the steady-state tolerances. ::
            %
            %     >>d.setSteadyTolerances(rtol, atol, component)
            %
            % :param rtol:
            %     Relative tolerance.
            % :param atol:
            %     Absolute tolerance.
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.

            if nargin < 4 | strcmp(component, 'default')
                ct.impl.call('mDomain_setSteadyTolerances', obj.domainID, rtol, atol, -1);
            elseif iscell(component)
                for ii = 1:length(component)
                    n = obj.componentIndex(component{ii});
                    ct.impl.call('mDomain_setSteadyTolerances', obj.domainID, rtol, atol, n);
                end
            else
                n = obj.componentIndex(component);
                ct.impl.call('mDomain_setSteadyTolerances', obj.domainID, rtol, atol, n);
            end

        end

        function setTransientTolerances(obj, rtol, atol, component)
            % Set the transient tolerances. ::
            %
            %     >> d.setTransientTolerances(rtol, atol, component)
            %
            % :param rtol:
            %     Relative tolerance.
            % :param atol:
            %     Absolute tolerance.
            % :param component:
            %     String or cell array of strings of component values
            %     whose tolerances should be set. If ``'default'`` is
            %     specified, the tolerance of all components will be set.

            if nargin < 4 | strcmp(component, 'default')
                ct.impl.call('mDomain_setTransientTolerances', obj.domainID, rtol, atol, -1);
            elseif iscell(component)
                for ii = 1:length(component)
                    n = obj.componentIndex(component{ii});
                    ct.impl.call('mDomain_setTransientTolerances', obj.domainID, rtol, atol, n);
                end
            else
                n = obj.componentIndex(component);
                ct.impl.call('mDomain_setTransientTolerances', obj.domainID, rtol, atol, n);
            end

        end

    end

end
