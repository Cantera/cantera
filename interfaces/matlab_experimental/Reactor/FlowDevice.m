classdef FlowDevice < handle

    properties
        type
        id
        upstream
        downstream
    end

    methods
        %% FlowDevice class constructor

        function x = FlowDevice(typ)
            % Flow Device class constructor.
            %
            % :parameter typ:
            %    Type of flow device to be created. Type =
            %    'MassFlowController', 'PressureController' or 'Valve'.

            checklib;

            if nargin == 0
                error('please specify the type of flow device to be created');
            end

            x.type = typ;
            x.id = calllib(ct, 'flowdev_new', typ);
%             if x.id < 0
%                 error(geterr);
%             end
            x.upstream = -1;
            x.downstream = -1;
        end

        %% Utility methods

        function clear(f)
            % Clear the specified flow device from memory.

            calllib(ct, 'flowdev_del', f.id);
        end

        %% FlowDevice methods

        function install(f, upstream, downstream)
            % Install a flow device between reactors or reservoirs.
            %
            % :parameter upstream:
            %    Upsteram 'Reactor' or 'Reservoir'.
            % :parameter downstream:
            %    Downstream 'Reactor' or 'Reservoir'.

            if nargin == 3
                if ~isa(upstream, 'Reactor') || ~isa(downstream, 'Reactor')
                    error(['Flow devices can only be installed between',...
                           'reactors or reservoirs']);
                end
                i = upstream.id;
                j = downstream.id;
                calllib(ct, 'flowdev_install', f.id, i, j);
%                 if ok < 0
%                     error(geterr)
%                 end
            else error('install requires 3 arguments');
            end
        end

        function mdot = massFlowRate(f)
            % Get the mass flow rate.
            %
            % :return:
            %    The mass flow rate through the flow device

            mdot = calllib(ct, 'flowdev_massFlowRate2', f.id);
        end

        function setFunction(f, mf)
            % Set the time function with class 'func'.
            %
            % :parameter mf:
            %    Instance of class 'func'.

            if strcmp(f.type, 'MassFlowController')
                k = calllib(ct, 'flowdev_setTimeFunction', f.id, ...
                            mf.id);
%                 if k < 0
%                     error(geterr);
%                 end
            else
                error('Time function can only be set for mass flow controllers.');
            end
        end

        function setMassFlowRate(f, mdot)
            % Set the mass flow rate to a constant value.
            %
            % :parameter mdot:
            %    Mass flow rate

            if strcmp(f.type, 'MassFlowController')
                k = calllib(ct, 'flowdev_setMassFlowCoeff', f.id, mdot);
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
            % is negative, zero is :returned.
            %
            % :parameter k:
            %    Value fo the valve coefficient. Unit: kg/Pa-s.

            if ~strcmp(f.type, 'Valve')
                error('Valve coefficient can only be set for valves.');
            end
            ok = calllib(ct, 'flowdev_setValveCoeff', f.id, k);
%            if k < 0
%                error(geterr);
%            end
        end

    end
end
