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

    properties (SetAccess = protected)

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

        function r = ReactorNet(reactors)
            % Create a :mat:class:`ReactorNet` object.

            ctIsLoaded;

            if nargin ~= 1
                error('Wrong number of arguments to ReactorNet constructor.');
            end

            if isa(reactors, 'Reactor')
                % Allow simpler syntax for creating a network with one
                % reactor.
                reactors = {reactor};
            end

            r.id = ctFunc('reactornet_new');

            % add reactors
            nr = length(reactors);

            for i = 1:nr
                ctFunc('reactornet_addreactor', r.id, reactors{i}.id);
            end

        end

        %% ReactorNet Class Destructor

        function delete(r)
            % Delete the :mat:class:`ReactorNet` object object.

            if isempty(r.id)
                return
            end
            ctFunc('reactornet_del', r.id);
        end

        %% ReactorNet Utility Methods

        function advance(r, tout)
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

            ctFunc('reactornet_advance', r.id, tout);
        end

        %% ReactorNet set methods

        function set.time(r, t)
            ctFunc('reactornet_setInitialTime', r.id, t);
        end

        function set.maxTimeStep(r, maxstep)
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

            ctFunc('reactornet_setMaxTimeStep', r.id, maxstep);
        end

        function setSensitivityTolerances(r, rerr, aerr)
            % Set the error tolerance for sensitivity analysis. ::
            %
            %     >> r.setSensitivityTOlerances(rerr, aerr)
            %
            % :param rerr:
            %    Scalar relative error tolerance.
            % :param aerr:
            %    Scalar absolute error tolerance.

            ctFunc('reactornet_setSensitivityTolerances', r.id, rerr, aerr);
        end

        function set.atol(r, aeer)
            rerr = r.rtol;
            ctFunc('reactornet_setTolerances', r.id, rerr, aerr);
        end

        function set.rtol(r, reer)
            aerr = r.atol;
            ctFunc('reactornet_setTolerances', r.id, rerr, aerr);
        end

        %% ReactorNet get methods

        function t = get.dt(r)
            t = ctFunc('reactor_step', r.id);
        end

        function t = get.time(r)
            t = ctFunc('reactornet_time', r.id);
        end

        function t = get.atol(r)
            t = ctFunc('reactornet_atol', r.id);
        end

        function t = get.rtol(r)
            t = ctFunc('reactornet_rtol', r.id);
        end

        function s = sensitivity(r, component, p, rxtr)
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
                ctFunc('reactornet_sensitivity', r.id, component, p, rxtr.id);
            end

        end

    end

end
