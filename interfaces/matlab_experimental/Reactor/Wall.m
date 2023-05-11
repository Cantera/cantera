classdef Wall < handle
    % Wall Class ::
    %
    %     >> x = Wall(l, r)
    %
    % A Wall separates two reactors, or a reactor and a reservoir.
    % A Wall has a finite area, may conduct heat between the two
    % reactors on either side, and may move like a piston.
    %
    % Walls are stateless objects in Cantera, meaning that no
    % differential equation is integrated to determine any wall
    % property. Since it is the wall (piston) velocity that enters
    % the energy equation, this means that it is the velocity, not
    % the acceleration or displacement, that is specified. The wall
    % velocity is computed from:
    %
    % .. math:: v = K(P_{left} - P_{right}) + v_0(t),
    %
    % where K is a non-negative constant, and v_0 is a specified
    % function of time. The velocity is positive if the wall is
    % moving to the right.
    %
    % The heat flux through the wall is computed from:
    %
    %  .. math:: q = U(T_{left} - T_{right}) + q_0(t),
    %
    % where U is the overall heat transfer coefficient for
    % conduction/convection, and q_0 is a specified function of
    % time. The heat flux is positive when heat flows from the
    % reactor on the left to the reactor on the right.
    %
    % Note: the Wall class constructor only assign default values
    % to various properties. The user could specify those properties
    % after initial construction by using the various methods of
    % the Wall class.
    %
    % :param l:
    %    Instance of class :mat:class:`Reactor` to be used as the bulk phase
    %    on the left side of the wall.
    % :param r:
    %    Instance of class :mat:class:`Reactor` to be used as the bulk phase
    %    on the right side of the wall.
    % :return:
    %    Instance of class :mat:class:`Wall`.

    properties (SetAccess = immutable)

        id
        type

    end

    properties (SetAccess = protected)

        left % Reactor on the left.
        right % Reactor on the right.

    end

    properties (SetAccess = public)

        area % Area of the wall in m^2.
        thermalResistance % Thermal resistance in K*m^2/W.
        heatTransferCoeff % Heat transfer coefficient in W/(m^2-K).
        emissivity % Non-dimensional emissivity.
        expansionRateCoeff % Expansion rate coefficient in m/(s-Pa).

        % Heat flux in W/m^2.
        %
        % Must be set by an instance of :mat:class:`Func`, which allows the
        % heat flux to be an arbitrary function of time. It is possible
        % to specify a constant heat flux by using the polynomial
        % functor with only the first term specified.
        heatFlux

        % Velocity in m/s.
        %
        % Must be set by an instance of :mat:class:`Func`, which allows the
        % velocity to be an arbitrary function of time. It is possible
        % to specify a constant velocity by using the polynomial
        % functor with only the first term specified.
        velocity
    end

    methods
        %% Wall Class Constructor

        function w = Wall(l, r)
            % Create a :mat:class:`Wall` object.
            ctIsLoaded;

            % At the moment, only one wall type is implemented
            typ = 'Wall';

            w.type = char(typ);
            w.id = ctFunc('wall_new', w.type);

            % Install the wall between left and right reactors
            w.left = l;
            w.right = r;
            ctFunc('wall_install', w.id, l.id, r.id);

            % Set default values.
            w.area = 1.0;
            w.expansionRateCoeff = 0.0;
            w.heatTransferCoeff = 0.0;

            % Check whether the wall is ready.
            ok = ctFunc('wall_ready', w.id);
            if ok
                disp('The wall object is ready.');
            else
                error('The wall object is not ready.');
            end

        end

        %% Wall Class Destructor

        function delete(w)
            % Clear the :mat:class:`Wall` object.

            ctFunc('wall_del', w.id);
        end

        %% ReactorNet get methods

        function a = get.area(w)
            a = ctFunc('wall_area', w.id);
        end

        function q = qdot(w, t)
            % Total heat transfer rate through a wall at a given time t.

            q = ctFunc('wall_Q', w.id, t);
        end

        function v = vdot(w, t)
            % Rate of volumetric change at a given time t.
            v = ctFunc('wall_vdot', w.id, t);
        end

        %% ReactorNet set methods

        function set.area(w, a)
            ctFunc('wall_setArea', w.id, a);
        end

        function set.thermalResistance(w, r)
            ctFunc('wall_setThermalResistance', w.id, r);
        end

        function set.heatTransferCoeff(w, u)
            ctFunc('wall_setHeatTransferCoeff', w.id, u);
        end

        function set.emissivity(w, epsilon)
            ctFunc('wall_setEmissivity', w.id, epsilon);
        end

        function set.expansionRateCoeff(w, k)
            ctFunc('wall_setExpansionRateCoeff', w.id, k);
        end

        function set.heatFlux(w, f)
            ctFunc('wall_setHeatFlux', w.id, f.id);
        end

        function set.velocity(w, f)
            ctFunc('wall_setVelocity', w.id, f.id);
        end

    end

end
