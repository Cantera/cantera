classdef FlowDevice < handle

    properties
        type
        id
        upstream
        downstream
    end

    properties(Constant = true)
        lib = 'cantera_shared'
    end

    methods
        %% FlowDevice class constructor

        function x = FlowDevice(typ)
            % Flow Device class constructor.
            %
            % :param typ:
            %    Type of flow device to be created. Type =
            %    'MassFlowController', 'PressureController' or 'Valve'.

            checklib;

            if nargin == 0
                typ = 'MassFlowController';
            end

            if isa(typ, 'double')
                warning(['Definition via integer type to be deprecated', ...
                         ' after Cantera 2.5.']);
                device_types = {'MassFlowController', 'PressureController', ...
                                'Valve'};
                typ = device_types(typ);
            end

            x.type = typ;
            x.id = calllib(x.lib, 'flowdev_new2', typ);
%             if x.id < 0
%                 error(geterr);
%             end
            x.upstream = -1;
            x.downstream = -1;
        end

        %% Utility methods

        function clear(f)
            % Clear the specified flow device from memory.
            checklib;
            calllib(f.lib, 'flowdev_del', f.id);
        end

        %% FlowDevice methods

        function install(f, upstream, downstream)
            % Install a flow device between reactors or reservoirs.
            %
            % :param upstream:
            %    Upsteram 'Reactor' or 'Reservoir'.
            % :param downstream:
            %    Downstream 'Reactor' or 'Reservoir'.

            checklib;
            if nargin == 3
                if ~isa(upstream, 'Reactor') || ~isa(downstream, 'Reactor')
                    error(['Flow devices can only be installed between',...
                           'reactors or reservoirs']);
                end
                i = upstream.id;
                j = downstream.id;
                ok = calllib(f.lib, 'flowdev_install', f.id, i, j);
%                 if ok < 0
%                     error(geterr)
%                 end
            else error('install requires 3 arguments');
            end
        end

        function mdot = massFlowRate(f, time)
            % Get the mass flow rate at a given time.
            %
            % :param time:
            %    Time at which the mass flow rate is desired.
            % :return:
            %    The mass flow rate through the flow device at the given
            %    time.

            checklib;
            if nargin == 1
                mdot = calllib(f.lib, 'flowdev_massFlowRate2', f.id);
            else
                warning(['"Time" argument to massFlowRate is deprecated', ...
                         'and will be removed after Cantera 2.5.']);
                mdot = calllib(f.lib, 'flowdev_massFlowRate', f.id, time);
            end
        end

        function setFunction(f, mf)
            % Set the time function with class 'func'.
            %
            % :param mf:
            %    Instance of class 'func'.

            checklib;
            if strcmp(f.type, 'MassFlowController')
                k = calllib(f.lib, 'flowdev_setTimeFunction', f.id, ...
                            mf.id);
%                 if k < 0
%                     error(geterr);
%                 end
            else
                error('Mass flow rate can only be set for mass flow controllers.');
            end
        end

        function setMassFlowRate(f, mdot)
            % Set the mass flow rate to a constant value.
            %
            % :param mdot:
            %    Mass flow rate

            checklib;
            if strcmp(f.type, 'MassFlowController')
                k = calllib(f.lib, 'flowdev_setMassFlowCoeff', f.id, mdot);
%                 if k < 0
%                     error(geterr);
%                 end
            else
                error('Mass flow rate can only be set for mass flow controllers.');
            end
        end

        function setValveCoeff(f, k)
            % Set the valve coefficient 'K'.
            %
            % The mass flow rate [kg/s] is computed from the expression
            % mdot = K(P_upstream - P_downstream)
            % as long as this produces a positive value. If thsi expression
            % is negative, zero is returned.
            %
            % :param k:
            %    Value fo the valve coefficient. Unit: kg/Pa-s.

            checklib;
            if ~strcmp(f.type, 'Valve')
                error('Valve coefficient can only be set for valves.');
            end
            ok = calllib(f.lib, 'flowdev_setValveCoeff', f.id, k);
%            if k < 0
%                error(geterr);
%            end
        end

    end
end
