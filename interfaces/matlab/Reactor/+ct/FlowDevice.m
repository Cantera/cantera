classdef (Abstract) FlowDevice < ct.Connector
    % FlowDevice Class ::
    %
    %     >> x = ct.FlowDevice(typ, name)
    %
    % Base class for devices that allow flow between reactors.
    % :mat:class:`ct.FlowDevice` objects are assumed to be adiabatic,
    % non-reactive, and have negligible internal volume, so that they are
    % internally always in steady-state even if the upstream and downstream
    % reactors are not. The fluid enthalpy, chemical composition, and mass
    % flow rate are constant across a :mat:class:`ct.FlowDevice`, and the
    % pressure difference equals the difference in pressure between the
    % upstream and downstream reactors.
    %
    % See also: :mat:class:`ct.MassFlowController`, :mat:class:`ct.Valve`
    %
    % :param typ:
    %     Type of :mat:class:`ct.FlowDevice` to be created. ``typ='MassFlowController'``
    %     for :mat:class:`ct.MassFlowController`,  ``typ='PressureController'`` for
    %     :mat:class:`ct.PressureController`, and ``typ='Valve'`` for
    %     :mat:class:`ct.Valve`.
    % :param upstream:
    %     Upstream reactor or reservoir.
    % :param downstream:
    %     Downstream Reactor or reservoir.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).

    properties (SetAccess = immutable)

        % Upstream object of type :mat:class:`ct.ReactorBase`.
        upstream

        % Downstream object of type :mat:class:`ct.ReactorBase`.
        downstream

    end

    properties (SetAccess = public)

        % The mass flow rate through the :mat:class:`ct.FlowDevice` at the current time.
        %
        % The setter method can either take a double value or a function represented by
        % an instance of :mat:class:`ct.Func1`.
        massFlowRate

        % Valve coefficient in kg/s/Pa.
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

        function obj = FlowDevice(typ, upstream, downstream, name)
            arguments
                typ (1,1) string
                upstream (1,1) ct.ReactorBase
                downstream (1,1) ct.ReactorBase
                name (1,1) string = "(none)"
            end

            obj@ct.Connector(typ, upstream, downstream, name);
            obj.upstream = upstream;
            obj.downstream = downstream;
        end

        %% FlowDevice Get Methods

        function mdot = get.massFlowRate(obj)
            mdot = ct.impl.call('mFlowdev_massFlowRate', obj.id);
        end

        %% FlowDevice Set Methods

        function set.massFlowRate(obj, mdot)

            if strcmp(obj.type, 'MassFlowController')
                if isa(mdot, 'double')
                    ct.impl.call('mFlowdev_setDeviceCoefficient', obj.id, mdot);
                elseif isa(mdot, 'ct.Func1')
                    ct.impl.call('mFlowdev_setTimeFunction', obj.id, mdot.id);
                else
                    error('Mass flow rate must either be a value or function.');
                end
            else
                error('Mass flow rate can only be set for mass flow controllers.');
            end

        end

        function setPrimary(obj, d)
            % Set the primary flow device used to compute this device's
            % mass flow rate. ::
            %
            %     >> f.setPrimary(d)
            %
            % :param d:
            %     Instance of class :mat:class:`ct.FlowDevice`.

            if strcmp(obj.type, 'PressureController')
                ct.impl.call('mFlowdev_setPrimary', obj.id, d);
            else
                error('Primary flow device can only be set for pressure controllers.');
            end

        end

        function set.valveCoeff(obj, k)

            if ~strcmp(obj.type, 'Valve')
                error('Valve coefficient can only be set for valves.');
            end

            ct.impl.call('mFlowdev_setDeviceCoefficient', obj.id, k);
        end

    end

end
