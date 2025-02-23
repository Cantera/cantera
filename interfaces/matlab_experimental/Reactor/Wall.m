classdef Wall < handle
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
    % :param name:
    %     Wall name (optional; default is ``(none)``).
    % :return:
    %    Instance of class :mat:class:`Wall`.

    properties (SetAccess = immutable)

        id
        type

    end

    properties (SetAccess = protected)

        name  % Name of wall.

        left % Reactor on the left.
        right % Reactor on the right.

    end

    properties (SetAccess = public)

        area % Area of the wall in m^2.
        heatRate % Total heat transfer rate through the wall at current time step in W.
        expansionRate % Rate of volumetric change at current time step in m^3/s.
        thermalResistance % Thermal resistance in K*m^2/W.
        heatTransferCoeff % Heat transfer coefficient in W/(m^2-K).
        emissivity % Non-dimensional emissivity.
        expansionRateCoeff % Expansion rate coefficient in m/(s-Pa).

        % Heat flux in W/m^2.
        %
        % Must be set by an instance of :mat:class:`Func1`, which allows the
        % heat flux to be an arbitrary function of time. It is possible
        % to specify a constant heat flux by using the polynomial
        % functor with only the first term specified.
        heatFlux

        % Velocity in m/s.
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
            % Create a :mat:class:`Wall` object.
            ctIsLoaded;

            % At the moment, only one wall type is implemented
            typ = 'Wall';
            if nargin < 3
                name = '(none)';
            end

            w.type = char(typ);
            w.id = ctFunc('wall_new', w.type, name);

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
            if ~ok
                error('The wall object is not ready.');
            end

        end

        %% Wall Class Destructor

        function delete(w)
            % Clear the :mat:class:`Wall` object.

            if isempty(w.id)
                return
            end
            ctFunc('wall_del', w.id);
        end

        %% ReactorNet get methods

        function name = get.name(w)
            name = ctString('wall_name', w.id);
        end

        function a = get.area(w)
            a = ctFunc('wall_area', w.id);
        end

        function q = get.heatRate(w)
            q = ctFunc('wall_heatRate', w.id);
        end

        function v = get.expansionRate(w)
            v = ctFunc('wall_expansionRate', w.id);
        end

        %% ReactorNet set methods

        function set.name(w, name)
            ctFunc('wall_setName', w.id, name);
        end

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
