classdef ReactorNet < handle

    properties
        id
        time
        dt
        atol
        rtol
    end

    methods
        %% ReactorNet class constructor

        function r = ReactorNet(reactors)
            % A 'ReactorNet' object is a container that holds one ormore
            % 'Reactor' objects in a network. 'ReactorNet' objects are used
            % to simultaneously advance the state of one or more coupled
            % reactors.
            %
            % parameter reactors:
            %    Instance of class 'Reactor' or a cell array of instance of
            %    'Reactor'.
            % return:
            %    Instance of class 'ReactorNet'.

            checklib;

            if nargin ~= 1
                error('Wrong number of arguments to ReactorNet constructor.');
            end

            if isa(reactors, 'Reactor')
                % Allow simpler syntax for creating a network with one
                % reactor.
                reactors = {reactor};
            end

            r.id = calllib(ct, 'reactornet_new');
%             if r.id < 0
%                 error(geterr);
%             end

            % add reactors
            nr = length(reactors);
            for i = 1:nr
                r.addReactor(reactors{i});
            end
        end

        %% Utility methods

        function clear(r)
            % Clear the ReactorNet object from the memory.
            checklib;
            calllib(ct, 'reactornet_del', r.id);
        end

        function addReactor(net, reactor)
            % Add a reactor to a network.
            %
            % parameter net:
            %    Instance of class 'ReactorNet'.
            % parameter reactor:
            %    Instance of class 'Solution'.

            checklib;
            calllib(ct, 'reactornet_addreactor', net.id, reactor.id);
        end

        function advance(r, tout)
            % Advance the state of the reactor network in time.
            %
            % Method 'advance' integrates the system of ODEs that determine
            % the rate of change of the volume, the mass of each species,
            % and the total energy for each reactor. The integration is
            % carried out from the current time to "tout". (Note: "tout" is
            % an absolute time, not a time interval.) The integrator may
            % take many internal timesteps before reaching tout.
            %
            % parameter tout:
            %    End time of the integration. Unit: s.

            checklib;
            calllib(ct, 'reactornet_advance', r.id, tout);
        end

        %% ReactorNet set methods

        function setInitialTime(r, t)
            % Set the initial time of the integration.
            %
            % parameter t:
            %    Time at which integration should be restarted, using the
            %    current state as the initial condition. Unit: s.

            checklib;
            calllib(ct, 'reactornet_setInitialTime', r.id, t);
        end

        function setMaxTimeStep(r, maxstep)
            % Set the maximum time step.
            %
            % The integrator chooses a step size based on the desired error
            % tolerance and the rate at which the solution is changing.
            % In some cases, the solution changes very slowly at first,
            % then very rapidly (ifnition problems). In such cases, the
            % integrator may choose a timestep that is too large, which
            % leads to numerical problems later. Use thismethod to set an
            % upper bound on the timestep.

            checklib;
            calllib(ct, 'reactornet_setMaxTimeStep', r.id, maxstep);
        end

        function setTolerances(r, rerr, aerr)
            % Set the error tolerance.
            %
            % parameter rtol:
            %    Scalar relative error tolerance.
            % parameter atol:
            %    Scalar absolute error tolerance.

            checklib;
            calllib(ct, 'reactornet_setTolerances', r.id, rerr, aerr);
        end

        %% ReactorNet get methods

        function t = get.dt(r)
            % Get the internal time step in s.
            %
            % The integrator used to integrate the ODEs (CVODE) takes
            % variable-size steps, chosen so that a specified error
            % tolerance is maintained. At times when the solution is
            % rapidly changing, the time step becomes smaller to resolve
            % the solution.
            checklib;
            t = calllib(ct, 'reactor_step', r.id);
        end

        function t = get.time(r)
            % Get the current time in s.
            checklib;
            t = calllib(ct, 'reactornet_time', r.id);
        end

        function t = get.rtol(r)
            % Get the relative error tolerance
            checklib;
            t = calllib(ct, 'reactornet_rtol', r.id);
        end

        function t = get.atol(r)
            % Get the absolute error tolerance
            checklib;
            t = calllib(ct, 'reactornet_atol', r.id);
        end

    end
end
