classdef Wall < handle

    properties
        id
        type
        left
        right
        area
    end

    methods
        %% Wall class constructor

        function x = Wall(left, right, area, k, u, q, v)
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
            % Note: all of the arguments are optional and can be activated
            % after initial construction by using the various methods of
            % the 'Wall' class. Any improperly specified arguments will
            % generate warnings; these can be ignored if the intention was
            % to not use a particular argument. Thus, the velocity of the
            % wall can be set by using empty strings or 0.0 for each of the
            % arguments before the velocity with no harm.
            %
            % :parameter left:
            %    Instance of class 'Reactor' to be used as the bulk phase
            %    on the left side of the wall.
            % :parameter right:
            %    Instance of class 'Reactor' to be used as the bulk phase
            %    on the right side of the wall.
            % :parameter area:
            %    The area of the wall in m^2. Defaults to 1.0 m^2 if not
            %    specified.
            % :parameter k:
            %    Expansion rate coefficient in m/(s-Pa). Defaults to 0.0 if
            %    not specified.
            % :parameter u:
            %    Heat transfer coefficient in W/(m^2-K). Defaults to 0.0 if
            %    not specified.
            % :parameter q:
            %    Heat flux in W/m^2. Defaults to 0.0 if not specified. Must
            %    be an instance of 'Func'.
            % :parameter v:
            %    Velocity of the wall in m/s. Defaults to 0.0 if not
            %    specified. Must be an instance of 'Func'.

            checklib;

            % At the moment, only one wall type is implemented
            typ = 'Wall';

            x.type = char(typ);
            x.id = callct('wall_new', x.type);
            if x.index < 0
                error(geterr);
            end
            x.left = -1;
            x.right = -1;

            if nargin >= 2
                if isa(left, 'Reactor') && isa(right, 'Reactor')
                    x.install(left, right);
                else
                    warning(['left and/or right were not instances of Reactor, ' ...
                             'and were not installed.'])
                end
            end

            if nargin >= 3
                if isnumeric(area)
                    x.setArea(area);
                else
                    warning('area was not a number and the area was not set')
                end
            end

            if nargin >= 4
                if isnumeric(k)
                    x.setExpansionRateCoeff(k);
                else
                    warning(['k was not a number and the expansion rate ',...
                             'coefficient was not set'])
                end
            end

            if nargin >= 5
                if isnumeric(u)
                    x.setHeatTransferCoeff(u);
                else
                    warning(['u was not a number and the expansion rate ',...
                             'coefficient was not set'])
                end
            end

            if nargin >= 6
                if isa(q, 'Func')
                    x.setHeatFlux(q);
                else
                    warning('q was not an instance of Func and was not set')
                end
            end

            if nargin >= 7
                if isa(v, 'Func')
                    x.setVelocity(v)
                else
                    warning('v was not an instance of Func and was not set')
                end
            end
        end

        %% Utility methods

        function clear(w)
            % Clear the Wall object from the memory.

            callct('wall_del', w.id);
        end

        function install(w, l, r)
            % Install a wall between two reactors.

            w.left = l;
            w.right = r;
            callct('wall_install', w.id, l.id, r.id);
        end

        function ok = ready(w)
            % Check whether a wall is ready.

            ok = callct('wall_ready', w.id);
        end

        %% ReactorNet set methods

        function set.area(w, a)
            % Set the area of a wall.

            callct('wall_setArea', w.id, a);
        end

        function setThermalResistance(w, r)
            % Set the thermal resistance.
            %
            % :parameter r:
            %    Thermal resistance. Unit: K*m^2/W.

            callct('wall_setThermalResistance', w.id, r);
        end

        function setHeatTransferCoeff(w, u)
            % Set the thermal transfer coefficient.
            %
            % :parameter u:
            %    Heat transfer coefficient. Unit: W/(m^2-K).

            callct('wall_setHeatTransferCoeff', w.id, u);
        end

        function setEmissivity(w, epsilon)
            % Set the emissivity.
            %
            % :param epsilon:
            %    Nondimensional emissivity.

            callct('wall_setEmissivity', w.id, epsilon);
        end

        function setExpansionRateCoeff(w, k)
            % Set the expansion rate coefficient.
            %
            % :parameter k:
            %    Expanstion rate coefficient. Unit: m/(s-Pa).

            callct('wall_setExpansionRateCoeff', w.id, k);
        end

        function setHeatFlux(w, f)
            % Set the heat flux.
            %
            % Must be set by an instance of class 'Func', which allows the
            % heat flux to be an arbitrary function of time. It is possible
            % to specify a constant heat flux by using the polynomial
            % functor with only the first term specified.
            %
            % :parameter f:
            %    Instance of class 'Func'. Unit: W/m^2.

            callct('wall_setHeatFlux', w.id, f.id);
        end

        function setVelocity(w, f)
            % Set the velocity of the wall.
            %
            % Must be set by an instance of class 'Func', which allows the
            % velocity to be an arbitrary function of time. It is possible
            % to specify a constant velocity by using the polynomial
            % functor with only the first term specified.
            %
            % :parameter f:
            %    Instance of class 'Func'. Unit: m/s.

            callct('wall_setVelocity', w.id, f.id);
        end

        %% ReactorNet get methods

        function a = get.area(w)
            % Get the area of the wall in m^2.

            a = callct('wall_area', w.id);
        end

        function q = qdot(w, t)
            % Get the total heat transfer through a wall at given time.

            q = callct('wall_Q', w.id, t);
        end

        function v = vdot(w, t)
            % Get the rate of volumetric change at a given time.

            v = callct('wall_vdot', w.id, t);
        end

    end
end
