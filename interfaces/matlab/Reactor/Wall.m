classdef Wall < Connector
    % Wall Class ::
    %
    %     >> x = Wall(l, r, name)
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
    % where K is a non-negative constant, and :math:`v_0` is a specified
    % function of time. The velocity is positive if the wall is
    % moving to the right.
    %
    % The heat flux through the wall is computed from:
    %
    % .. math:: q = U(T_{left} - T_{right}) + q_0(t),
    %
    % where U is the overall heat transfer coefficient for
    % conduction/convection, and :math:`q_0` is a specified function of
    % time. The heat flux is positive when heat flows from the
    % reactor on the left to the reactor on the right.
    %
    % Note: the Wall class constructor only assign default values
    % to various properties. The user could specify those properties
    % after initial construction by using the various methods of
    % the Wall class.
    %
    % :param l:
    %    Instance of class :mat:class:`ReactorBase` to be used as the bulk phase
    %    on the left side of the wall.
    % :param r:
    %    Instance of class :mat:class:`ReactorBase` to be used as the bulk phase
    %    on the right side of the wall.
    % :param name:
    %     Wall name (optional; default is ``(none)``).

    properties (SetAccess = immutable)

        left % Reactor on the left.
        right % Reactor on the right.

    end

    properties (SetAccess = public)

        area % Area of the wall [m²].
        heatRate % Total heat transfer rate [W] through the wall at the current time.
        expansionRate % Rate of volumetric change [m³/s] at the current time.
        thermalResistance % Thermal resistance [K·m²/W].
        heatTransferCoeff % Heat transfer coefficient [W/m²/K].
        emissivity % Emissivity [-].
        expansionRateCoeff % Expansion rate coefficient [m/s/Pa].

        % Heat flux [W/m²].
        %
        % Must be set by an instance of :mat:class:`Func1`, which allows the
        % heat flux to be an arbitrary function of time. It is possible
        % to specify a constant heat flux by using the polynomial
        % functor with only the first term specified.
        heatFlux

        % Velocity [m/s].
        %
        % Must be set by an instance of :mat:class:`Func1`, which allows the
        % velocity to be an arbitrary function of time. It is possible
        % to specify a constant velocity by using the polynomial
        % functor with only the first term specified.
        velocity
    end

    methods
        %% Wall Class Constructor

        function w = Wall(l, r, name)
            arguments
                l (1,1) ReactorBase
                r (1,1) ReactorBase
                name (1,1) string = "(none)"
            end

            % Install the wall between left and right reactors
            w@Connector('Wall', l, r, name);
            w.left = l;
            w.right = r;

            % Set default values.
            w.area = 1.0;
            w.expansionRateCoeff = 0.0;
            w.heatTransferCoeff = 0.0;

        end

        %% Wall get methods

        function a = get.area(obj)
            a = ctFunc('mWall_area', obj.id);
        end

        function q = get.heatRate(obj)
            q = ctFunc('mWall_heatRate', obj.id);
        end

        function v = get.expansionRate(obj)
            v = ctFunc('mWall_expansionRate', obj.id);
        end

        %% Wall set methods

        function set.area(obj, a)
            ctFunc('mWall_setArea', obj.id, a);
        end

        function set.thermalResistance(obj, r)
            ctFunc('mWall_setThermalResistance', obj.id, r);
        end

        function set.heatTransferCoeff(obj, u)
            ctFunc('mWall_setHeatTransferCoeff', obj.id, u);
        end

        function set.emissivity(obj, epsilon)
            ctFunc('mWall_setEmissivity', obj.id, epsilon);
        end

        function set.expansionRateCoeff(obj, k)
            ctFunc('mWall_setExpansionRateCoeff', obj.id, k);
        end

        function set.heatFlux(obj, f)
            ctFunc('mWall_setHeatFlux', obj.id, f.id);
        end

        function set.velocity(obj, f)
            ctFunc('mWall_setVelocity', obj.id, f.id);
        end

    end

end
