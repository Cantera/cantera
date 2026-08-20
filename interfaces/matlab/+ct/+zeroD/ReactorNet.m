classdef ReactorNet < handle
    % ReactorNet class ::
    %
    %     >> r = ct.zeroD.ReactorNet(reactors)
    %
    % A :mat:class:`ct.zeroD.ReactorNet` object is a container that holds one or more
    % :mat:class:`ct.zeroD.ReactorBase` objects in a network. :mat:class:`ct.zeroD.ReactorNet`
    % objects are used to simultaneously advance the state of one or more coupled
    % reactors.
    %
    % :param reactors:
    %    An instance of or a cell array of instances of class
    %    :mat:class:`ct.zeroD.ReactorBase`.

    properties (SetAccess = immutable)
        id = -1

        % Cell array of the :mat:class:`ct.zeroD.ReactorBase` objects that make up
        % the network, in the order in which they were installed.
        reactors
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

        function obj = ReactorNet(reactors)
            arguments
                reactors
            end
            if isa(reactors, 'ct.zeroD.ReactorBase')
                % Allow simpler syntax for creating a network with one reactor.
                reactors = {reactors};
            end
            reactorIDs = cellfun(@(r) r.id, reactors);

            obj.id = ct.impl.call('mReactornet_new', reactorIDs);
            obj.reactors = reactors;
            obj.time = 0;

        end

        %% ReactorNet Class Destructor

        function delete(obj)
            % Delete the :mat:class:`ct.zeroD.ReactorNet` object object.
            if obj.id >= 0
                ct.impl.call('mReactornet_del', obj.id);
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

            ct.impl.call('mReactornet_advance', obj.id, tout);
        end

        function t = step(obj)
            % Take a single internal step. ::
            %
            %     >> n.step()
            %
            % The time/distance after taking the step is returned.

            t = ct.impl.call('mReactornet_step', obj.id);
        end

        %% ReactorNet set methods

        function set.time(obj, t)
            ct.impl.call('mReactornet_setInitialTime', obj.id, t);
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

            ct.impl.call('mReactornet_setMaxTimeStep', obj.id, maxstep);
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

            ct.impl.call('mReactornet_setSensitivityTolerances', obj.id, rerr, aerr);
        end

        function set.atol(obj, aerr)
            ct.impl.call('mReactornet_setAbsoluteTolerance', obj.id, aerr);
        end

        function set.rtol(obj, rerr)
            ct.impl.call('mReactornet_setRelativeTolerance', obj.id, rerr);
        end

        %% ReactorNet get methods

        function t = get.time(obj)
            t = ct.impl.call('mReactornet_time', obj.id);
        end

        function t = get.atol(obj)
            t = ct.impl.call('mReactornet_atol', obj.id);
        end

        function t = get.rtol(obj)
            t = ct.impl.call('mReactornet_rtol', obj.id);
        end

        function s = sensitivity(obj, component, p, r)
            % Sensitivity of a solution variable with respect to a sensitivity
            % parameter. ::
            %
            %     >> s = n.sensitivity(component, p)
            %     >> s = n.sensitivity(component, p, r)
            %
            % :param component:
            %    Name of the solution variable, given as a character array or a
            %    string, for example ``'temperature'`` or a species name. It is
            %    resolved within the reactor `r`. Alternatively, an integer giving
            %    the 1-based index of the variable in the global state vector of
            %    the network, in which case `r` is not used.
            % :param p:
            %    Integer, 1-based index of the sensitivity parameter.
            % :param r:
            %    The reactor in which `component` is to be found, given either as an
            %    instance of :mat:class:`ct.zeroD.ReactorBase` belonging to this
            %    network or as its 1-based position within the network. Defaults to
            %    the first reactor in the network.
            % :return:
            %    Scalar normalized sensitivity coefficient.

            arguments
                obj
                component
                p (1,1) double
                r = 1
            end

            if isnumeric(component)
                % Global index into the network state vector
                s = ct.impl.call('mReactornet_sensitivity', obj.id, ...
                                 component - 1, p - 1);
                return
            end

            if ~(ischar(component) || isstring(component))
                error(['Sensitivity component must be a character array, a ' ...
                       'string, or an integer index, not a %s.'], class(component));
            end

            if isa(r, 'ct.zeroD.ReactorBase')
                index = find(cellfun(@(x) x == r, obj.reactors), 1);
                if isempty(index)
                    error('The given reactor is not part of this network.');
                end
            elseif isnumeric(r) && isscalar(r)
                index = r;
            else
                error(['Reactor must be a ct.zeroD.ReactorBase object or a ' ...
                       '1-based index, not a %s.'], class(r));
            end

            s = ct.impl.call('mReactornet_sensitivityByName', obj.id, ...
                             char(component), p - 1, index - 1);

        end

    end

end
