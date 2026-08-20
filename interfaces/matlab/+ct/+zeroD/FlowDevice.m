classdef (Abstract) FlowDevice < ct.zeroD.Connector
    % FlowDevice Class.
    %
    % Base class for devices that allow flow between reactors. The constructor is called
    % by derived classes and cannot be used directly.
    %
    % :mat:class:`ct.zeroD.FlowDevice` objects are assumed to be adiabatic,
    % non-reactive, and have negligible internal volume, so that they are
    % internally always in steady-state even if the upstream and downstream
    % reactors are not. The fluid enthalpy, chemical composition, and mass
    % flow rate are constant across a :mat:class:`ct.zeroD.FlowDevice`, and the
    % pressure difference equals the difference in pressure between the
    % upstream and downstream reactors.
    %
    % See also: :mat:class:`ct.zeroD.MassFlowController`,
    % :mat:class:`ct.zeroD.PressureController`, :mat:class:`ct.zeroD.Valve`
    %
    % :param typ:
    %     Type of :mat:class:`ct.zeroD.FlowDevice` to be created. ``typ='MassFlowController'``
    %     for :mat:class:`ct.zeroD.MassFlowController`,  ``typ='PressureController'`` for
    %     :mat:class:`ct.zeroD.PressureController`, and ``typ='Valve'`` for
    %     :mat:class:`ct.zeroD.Valve`.
    % :param upstream:
    %     Upstream reactor or reservoir.
    % :param downstream:
    %     Downstream Reactor or reservoir.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).

    properties (SetAccess = immutable)

        % Upstream object of type :mat:class:`ct.zeroD.ReactorBase`.
        upstream

        % Downstream object of type :mat:class:`ct.zeroD.ReactorBase`.
        downstream

    end

    properties (SetAccess = public)

        % The mass flow rate through the :mat:class:`ct.zeroD.FlowDevice` at the current time.
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

        % The coefficient defining the behavior of this
        % :mat:class:`ct.zeroD.FlowDevice`, with a meaning that depends on the type
        % of the device: the mass flow rate [kg/s] for a
        % :mat:class:`ct.zeroD.MassFlowController`, the valve coefficient [kg/s/Pa]
        % for a :mat:class:`ct.zeroD.Valve`, and the pressure coefficient [kg/s/Pa]
        % for a :mat:class:`ct.zeroD.PressureController`.
        deviceCoefficient
    end

    methods
        %% FlowDevice Class Constructor

        function obj = FlowDevice(typ, upstream, downstream, name)
            arguments
                typ (1,1) string
                upstream (1,1) ct.zeroD.ReactorBase
                downstream (1,1) ct.zeroD.ReactorBase
                name (1,1) string = "(none)"
            end

            obj@ct.zeroD.Connector(typ, upstream, downstream, name);
            obj.upstream = upstream;
            obj.downstream = downstream;
        end

        %% FlowDevice Get Methods

        function mdot = get.massFlowRate(obj)
            mdot = ct.impl.call('mFlowdev_massFlowRate', obj.id);
        end

        function c = get.deviceCoefficient(obj)
            c = ct.impl.call('mFlowdev_deviceCoefficient', obj.id);
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
            %     Instance of class :mat:class:`ct.zeroD.FlowDevice`.

            if ~isa(d, 'ct.zeroD.FlowDevice')
                error('Primary device must be an instance of ct.zeroD.FlowDevice.');
            end

            if strcmp(obj.type, 'PressureController')
                ct.impl.call('mFlowdev_setPrimary', obj.id, d.id);
            else
                error('Primary flow device can only be set for pressure controllers.');
            end

        end

        function set.valveCoeff(obj, k)

            if ~strcmp(obj.type, 'Valve')
                error('Valve coefficient can only be set for valves.');
            end

            obj.deviceCoefficient = k;
        end

        function set.deviceCoefficient(obj, c)
            ct.impl.call('mFlowdev_setDeviceCoefficient', obj.id, c);
        end

    end

end
