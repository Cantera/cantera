classdef Wall < handle
    % Wall Class
    %
    % x = Wall()
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
    %    v = K(P_left - P_right) + v_0(t),
    %
    % where K is a non-negative constant, and v_0 is a specified
    % function of time. The velocity is positive if the wall is
    % moving to the right.
    %
    % The heat flux through the wall is computed from:
    %
    %    q = U(T_left - T_right) + q_0(t),
    %
    % where U is the overall heat transfer coefficient for
    % conduction/convection, and q_0 is a specified function of
    % time. The heat flux is positive when heat flows from the
    % reactor on the left to the reactor on the right.
    %
    % Note: the Wall class constructor only assign default values
    % to various properties. The user could specify those properties
    % after initial construction by using the various methods of
    % the 'Wall' class.
    %

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
    end

    methods
        %% Wall Class Constructor

        function x = Wall()
            checklib;

            % At the moment, only one wall type is implemented
            typ = 'Wall';

            x.type = char(typ);
            x.id = callct('wall_new', x.type);

            % Set default values.
            x.left = -1;
            x.right = -1;
            x.area = 1.0;
            x.setExpansionRateCoeff(0.0);
            x.setHeatTransferCoeff(0.0);

        end

        %% Wall Class Destructor

        function delete(w)
            % Clear the Wall object from the memory.

            callct('wall_del', w.id);
        end

        %% Wall Class Utility Methods

        function install(w, l, r)
            % Install a wall between two reactors.
            %
            % w.install(l, r)
            %
            % :param l:
            %    Instance of class 'Reactor' to be used as the bulk phase
            %    on the left side of the wall.
            % :param r:
            %    Instance of class 'Reactor' to be used as the bulk phase
            %    on the right side of the wall.

            w.left = l;
            w.right = r;
            callct('wall_install', w.id, l.id, r.id);
        end

        function ok = ready(w)
            % Check whether a wall is ready.
            %
            % ok = w.ready
            %

            ok = callct('wall_ready', w.id);
        end

        %% ReactorNet get methods

        function a = get.area(w)
            a = callct('wall_area', w.id);
        end

        function q = qdot(w, t)
            % Total heat transfer through a wall at a given time t.

            q = callct('wall_Q', w.id, t);
        end

        function v = vdot(w, t)
            % Rate of volumetric change at a given time t.
            v = callct('wall_vdot', w.id, t);
        end

        %% ReactorNet set methods

        function set.area(w, a)
            callct('wall_setArea', w.id, a);
        end

        function setThermalResistance(w, r)
            % Set the thermal resistance.
            %
            % w.setTHermalResistance(w)
            %
            % :param r:
            %    Thermal resistance. Unit: K*m^2/W.
            %

            callct('wall_setThermalResistance', w.id, r);
        end

        function setHeatTransferCoeff(w, u)
            % Set the thermal transfer coefficient.
            %
            % w.setHeatTransferCoeff(u)
            %
            % :param u:
            %    Heat transfer coefficient. Unit: W/(m^2-K).
            %

            callct('wall_setHeatTransferCoeff', w.id, u);
        end

        function setEmissivity(w, epsilon)
            % Set the emissivity.
            %
            % w.setEmissitivity(epsilon)
            %
            % :param epsilon:
            %    Nondimensional emissivity.
            %

            callct('wall_setEmissivity', w.id, epsilon);
        end

        function setExpansionRateCoeff(w, k)
            % Set the expansion rate coefficient.
            %
            % w.setExpansionRateCoeff(k)
            %
            % :param k:
            %    Expanstion rate coefficient. Unit: m/(s-Pa).
            %

            callct('wall_setExpansionRateCoeff', w.id, k);
        end

        function setHeatFlux(w, f)
            % Set the heat flux.
            %
            % w.setHeatFlux(f)
            %
            % Must be set by an instance of class 'Func', which allows the
            % heat flux to be an arbitrary function of time. It is possible
            % to specify a constant heat flux by using the polynomial
            % functor with only the first term specified.
            %
            % :param f:
            %    Instance of class 'Func'. Unit: W/m^2.
            %

            callct('wall_setHeatFlux', w.id, f.id);
        end

        function setVelocity(w, f)
            % Set the velocity of the wall.
            %
            % w.setVelocity(f)
            %
            % Must be set by an instance of class 'Func', which allows the
            % velocity to be an arbitrary function of time. It is possible
            % to specify a constant velocity by using the polynomial
            % functor with only the first term specified.
            %
            % :param f:
            %    Instance of class 'Func'. Unit: m/s.
            %

            callct('wall_setVelocity', w.id, f.id);
        end

    end

end
