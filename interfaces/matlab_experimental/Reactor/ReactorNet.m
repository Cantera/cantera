classdef ReactorNet < handle
    % ReactorNet class ::
    %
    %     >> r = ReactorNet(reactors)
    %
    % A :mat:class:`ReactorNet` object is a container that holds one or more
    % :mat:class:`Reactor` objects in a network. :mat:class:`ReactorNet` objects
    % are used to simultaneously advance the state of one or more coupled reactors.
    %
    % :param reactors:
    %    An instance of or a cell array of instances of class :mat:class:`Reactor`.
    % :return:
    %    Instance of class :mat:class:`ReactorNet`.

    properties (SetAccess = immutable)
        id
    end

    properties (SetAccess = public)

        % Internal time step in s. ::
        %
        %     >> t = r.dt
        %
        % The integrator used to integrate the ODEs (CVODE) takes
        % variable-size steps, chosen so that a specified error
        % tolerance is maintained. At times when the solution is
        % rapidly changing, the time step becomes smaller to resolve
        % the solution.
        dt

        % Current time in s.
        % The setter method sets the time at which integration should be
        % restarted, using the current state as the initial condition.
        time

        % Max time step in s.
        %
        % The integrator chooses a step size based on the desired error
        % tolerance and the rate at which the solution is changing.
        % In some cases, the solution changes very slowly at first,
        % then very rapidly (ignition problems). In such cases, the
        % integrator may choose a timestep that is too large, which
        % leads to numerical problems later. Use this method to set an
        % upper bound on the timestep.
        maxTimeStep

        atol % Absolute error tolerance.

        rtol % Relative error tolerance.

    end

    methods
        %% ReactorNet Class Constructor

        function n = ReactorNet(r)
            % Create a :mat:class:`ReactorNet` object.

            ctIsLoaded;

            if nargin == 0
                r = {};
            end

            if isa(r, 'Reactor')
                % Allow simpler syntax for creating a network with one
                % reactor.
                r = {r};
            end

            n.id = ctFunc('reactornet_new');

            % add reactors
            nr = length(r);

            for i = 1:nr
                ctFunc('reactornet_addreactor', n.id, r{i}.id);
            end

        end

        %% ReactorNet Class Destructor

        function delete(n)
            % Delete the :mat:class:`ReactorNet` object object.

            if isempty(n.id)
                return
            end
            ctFunc('reactornet_del', n.id);
        end

        %% ReactorNet Utility Methods

        function addReactor(n, r)
            % Add a reactor to the network.

            ctFunc('reactornet_addreactor', n.id, r.id);
        end

        function advance(n, tout)
            % Advance the state of the reactor network in time. ::
            %
            %     >> r.advance(tout)
            %
            % Method `advance` integrates the system of ODEs that determine
            % the rate of change of the volume, the mass of each species,
            % and the total energy for each reactor. The integration is
            % carried out from the current time to "tout". (Note: "tout" is
            % an absolute time, not a time interval.) The integrator may
            % take many internal timesteps before reaching tout.
            %
            % :param tout:
            %    End time of the integration. Unit: s.

            ctFunc('reactornet_advance', n.id, tout);
        end

        %% ReactorNet set methods

        function set.time(n, t)
            ctFunc('reactornet_setInitialTime', n.id, t);
        end

        function set.maxTimeStep(n, maxstep)
            % Set the maximum time step. ::
            %
            %     >> r.setMaxTimeStep(maxstep)
            %
            % The integrator chooses a step size based on the desired error
            % tolerance and the rate at which the solution is changing.
            % In some cases, the solution changes very slowly at first,
            % then very rapidly (ignition problems). In such cases, the
            % integrator may choose a timestep that is too large, which
            % leads to numerical problems later. Use this method to set an
            % upper bound on the timestep.
            %
            % :param maxstep:
            %    Scalar max time step.

            ctFunc('reactornet_setMaxTimeStep', n.id, maxstep);
        end

        function setSensitivityTolerances(n, rerr, aerr)
            % Set the error tolerance for sensitivity analysis. ::
            %
            %     >> r.setSensitivityTOlerances(nerr, aerr)
            %
            % :param rerr:
            %    Scalar relative error tolerance.
            % :param aerr:
            %    Scalar absolute error tolerance.

            ctFunc('reactornet_setSensitivityTolerances', n.id, rerr, aerr);
        end

        function set.atol(n, aeer)
            rerr = r.rtol;
            ctFunc('reactornet_setTolerances', n.id, rerr, aerr);
        end

        function set.rtol(n, reer)
            aerr = r.atol;
            ctFunc('reactornet_setTolerances', n.id, rerr, aerr);
        end

        %% ReactorNet get methods

        function t = get.dt(n)
            t = ctFunc('reactornet_step', n.id);
        end

        function t = get.time(n)
            t = ctFunc('reactornet_time', n.id);
        end

        function t = get.atol(n)
            t = ctFunc('reactornet_atol', n.id);
        end

        function t = get.rtol(n)
            t = ctFunc('reactornet_rtol', n.id);
        end

        function s = sensitivity(n, component, p, r)
            % Sensitivity of the solution variable `c` in reactor `rxtr`
            % with respect to the parameter `p` ::
            %
            %     >> s = r.sensitivity(component, p, rxtr)
            %
            % :param component:
            %    String name of variable.
            % :param p:
            %    Integer sensitivity parameter.
            % :param rxtr:
            %    Instance of class :mat:class:`reactor`.

            if isa(component, 'string')
                s = ctFunc('reactornet_sensitivity', n.id, component, p, r.id);
            end

        end

    end

end
