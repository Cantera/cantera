classdef FlowDevice < handle
    % FlowDevice Class ::
    %
    %     >> x = FlowDevice(typ, name)
    %
    % Base class for devices that allow flow between reactors.
    % :mat:class:`FlowDevice` objects are assumed to be adiabatic,
    % non-reactive, and have negligible internal volume, so that they are
    % internally always in steady-state even if the upstream and downstream
    % reactors are not. The fluid enthalpy, chemical composition, and mass
    % flow rate are constant across a :mat:class:`FlowDevice`, and the
    % pressure difference equals the difference in pressure between the
    % upstream and downstream reactors.
    %
    % See also: :mat:class:`MassFlowController`, :mat:class:`Valve`
    %
    % :param typ:
    %     Type of :mat:class:`FlowDevice` to be created. ``typ='MassFlowController'``
    %     for :mat:class:`MassFlowController`,  ``typ='PressureController'`` for
    %     :mat:class:`PressureController`, and ``typ='Valve'`` for :mat:class:`Valve`.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`FlowDevice`.

    properties (SetAccess = immutable)

        type % Type of flow device.
        id % ID of FlowDevice object.

    end

    properties (SetAccess = public)

        name  % Name of flow device.

        % Upstream object of type :mat:class:`Reactor` or :mat:class:`Reservoir`.
        upstream

        % Downstream object of type :mat:class:`Reactor` or :mat:class:`Reservoir`.
        downstream

        % The mass flow rate through the :mat:class:`FlowDevice` at the current time.
        %
        % The setter method can either take a double value or a function represented by
        % an instance of :mat:class:`Func1`.
        massFlowRate

        % Valve coefficient in kg/Pa-s.
        %
        % The mass flow rate [kg/s] is computed from the expression:
        %
        % .. math:: \dot{m} = K(P_{upstream} - P_{downstream})
        %
        % as long as this produces a positive value.  If this expression is
        % negative, zero is returned.
        valveCoeff
    end

    methods
        %% FlowDevice Class Constructor

        function x = FlowDevice(typ, name)
            % Create a :mat:class:`FlowDevice` object.

            ctIsLoaded;

            if nargin == 0
                error('please specify the type of flow device to be created');
            end
            if nargin < 2
                name = '(none)';
            end

            x.type = typ;
            x.id = ctFunc('flowdev_new', typ, name);
            x.upstream = -1;
            x.downstream = -1;
        end

        %% FlowDevice Class Destructor

        function delete(f)
            % Delete the :mat:class:`FlowDevice` object.

            if isempty(f.id)
                return
            end
            ctFunc('flowdev_del', f.id);
        end

        %% Utility Methods

        function install(f, upstream, downstream)
            % Install a flow device between reactors or reservoirs. ::
            %
            %     >> f.install(upstream, downstream)
            %
            % :param f:
            %     Instance of class :mat:class:`FlowDevice` to install.
            % :param upstream:
            %     Upstream :mat:class:`Reactor` or :mat:class:`Reservoir`.
            % :param downstream:
            %     Downstream :mat:class:`Reactor` or :mat:class:`Reservoir`.
            % :return:
            %     Instance of class :mat:class:`FlowDevice`.

            if nargin == 3

                if ~isa(upstream, 'Reactor') || ~isa(downstream, 'Reactor')
                    error(['Flow devices can only be installed between', ...
                           'reactors or reservoirs']);
                end

                i = upstream.id;
                j = downstream.id;
                ctFunc('flowdev_install', f.id, i, j);
            else error('install requires 3 arguments');
            end

        end

        %% FlowDevice Get Methods

        function name = get.name(f)
            name = ctString('flowdev_name', f.id);
        end

        function mdot = get.massFlowRate(f)
            mdot = ctFunc('flowdev_massFlowRate2', f.id);
        end

        %% FlowDevice Set Methods

        function set.name(f, name)
            ctFunc('flowdev_setName', f.id, name);
        end

        function set.massFlowRate(f, mdot)

            if strcmp(f.type, 'MassFlowController')
                if isa(mdot, 'double')
                    k = ctFunc('flowdev_setMassFlowCoeff', f.id, mdot);
                elseif isa(mdot, 'Func1')
                    k = ctFunc('flowdev_setTimeFunction', f.id, mdot.id);
                else
                    error('Mass flow rate must either be a value or function.');
                end
            else
                error('Mass flow rate can only be set for mass flow controllers.');
            end

        end

        function setPrimary(f, d)
            % Set the primary flow device used to compute this device's
            % mass flow rate. ::
            %
            %     >> f.setPrimary(d)
            %
            % :param f:
            %     Instance of class :mat:class:`MassFlowController`.
            % :param mf:
            %     Instance of class :mat:class:`Func1`.

            if strcmp(f.type, 'PressureController')
                k = ctFunc('flowdev_setPrimary', f.id, d);
            else
                error('Primary flow device can only be set for pressure controllers.');
            end

        end

        function set.valveCoeff(f, k)

            if ~strcmp(f.type, 'Valve')
                error('Valve coefficient can only be set for valves.');
            end

            ok = ctFunc('flowdev_setValveCoeff', f.id, k);
        end

    end

end
