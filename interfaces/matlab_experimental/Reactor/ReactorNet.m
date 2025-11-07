classdef ReactorNet < handle
    % ReactorNet class ::
    %
    %     >> r = ReactorNet(reactors)
    %
    % A :mat:class:`ReactorNet` object is a container that holds one or more
    % :mat:class:`ReactorBase` objects in a network. :mat:class:`ReactorNet` objects
    % are used to simultaneously advance the state of one or more coupled reactors.
    %
    % :param reactors:
    %    An instance of or a cell array of instances of class :mat:class:`ReactorBase`.

    properties (SetAccess = immutable)
        id = -1
    end

    properties (SetAccess = public)

        % Current time [s].
        % The setter method sets the time at which integration should be
        % restarted, using the current state as the initial condition.
        time

        % Max time step [s].
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

        function n = ReactorNet(reactors)
            arguments
                reactors
            end
            if isa(reactors, 'ReactorBase')
                % Allow simpler syntax for creating a network with one reactor.
                reactors = {reactors};
            end
            reactorIDs = cellfun(@(r) r.id, reactors);

            n.id = ctFunc('mReactornet_new', reactorIDs);
            n.time = 0;

        end

        %% ReactorNet Class Destructor

        function delete(obj)
            % Delete the :mat:class:`ReactorNet` object object.
            if obj.id >= 0
                ctFunc('mReactornet_del', obj.id);
            end
        end

        %% ReactorNet Utility Methods

        function advance(obj, tout)
            % Advance the state of the reactor network in time. ::
            %
            %     >> n.advance(tout)
            %
            % Method `advance` integrates the system of ODEs that determine
            % the rate of change of the volume, the mass of each species,
            % and the total energy for each reactor. The integration is
            % carried out from the current time to "tout". (Note: "tout" is
            % an absolute time, not a time interval.) The integrator may
            % take many internal timesteps before reaching tout.
            %
            % :param tout:
            %    End time [s] of the integration.

            ctFunc('mReactornet_advance', obj.id, tout);
        end

        function t = step(obj)
            % Take a single internal step. ::
            %
            %     >> n.step()
            %
            % The time/distance after taking the step is returned.

            t = ctFunc('mReactornet_step', obj.id);
        end

        %% ReactorNet set methods

        function set.time(obj, t)
            ctFunc('mReactornet_setInitialTime', obj.id, t);
        end

        function set.maxTimeStep(obj, maxstep)
            % Set the maximum time step. ::
            %
            %     >> n.setMaxTimeStep(maxstep)
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
            %    max time step [s].

            ctFunc('mReactornet_setMaxTimeStep', obj.id, maxstep);
        end

        function setSensitivityTolerances(obj, rerr, aerr)
            % Set the error tolerance for sensitivity analysis. ::
            %
            %     >> n.setSensitivityTOlerances(nerr, aerr)
            %
            % :param rerr:
            %    Scalar relative error tolerance.
            % :param aerr:
            %    Scalar absolute error tolerance.

            ctFunc('mReactornet_setSensitivityTolerances', obj.id, rerr, aerr);
        end

        function set.atol(obj, aerr)
            rerr = obj.rtol;
            ctFunc('mReactornet_setTolerances', obj.id, rerr, aerr);
        end

        function set.rtol(obj, rerr)
            aerr = obj.atol;
            ctFunc('mReactornet_setTolerances', obj.id, rerr, aerr);
        end

        %% ReactorNet get methods

        function t = get.time(obj)
            t = ctFunc('mReactornet_time', obj.id);
        end

        function t = get.atol(obj)
            t = ctFunc('mReactornet_atol', obj.id);
        end

        function t = get.rtol(obj)
            t = ctFunc('mReactornet_rtol', obj.id);
        end

        function s = sensitivity(obj, component, p, r)
            % Sensitivity of the solution variable `c` in reactor `r`
            % with respect to the parameter `p` ::
            %
            %     >> s = n.sensitivity(component, p, r)
            %
            % :param component:
            %    String name of variable.
            % :param p:
            %    Integer sensitivity parameter.
            % :param r:
            %    Instance of class :mat:class:`ReactorBase`.

            if isa(component, 'string')
                s = ctFunc('mReactornet_sensitivity', obj.id, component, p, r.id);
            end

        end

    end

end
