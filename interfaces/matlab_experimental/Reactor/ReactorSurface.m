classdef ReactorSurface < handle

    properties
        surfID
        area
        reactor
    end

    methods
        %% ReactorSurface class constructor

        function s = ReactorSurface(kleft, reactor, area)
            % A surface on which heterogeneous reactions take place. The
            % mechanism object (typically an instance of class 'Interface')
            % mustb be constructed so that it is properly linked to the
            % object representing the fluid in the reactor. The surface
            % temperature on each side is taken to be equal to the
            % temprature of the reactor.
            %
            % Note: all of the arguments are optional and can be activated
            % after initial construction by using the various methods of
            % the 'ReactorSurface' class.
            %
            % :parameter kleft:
            %    Surface reaction mechanisms for the left-facing surface.
            %    This must bean instance of class 'Kinetics', or of a class
            %    derived from 'Kinetics', such as 'Interface'.
            % :parameter reactor:
            %    Instance of class 'Reactor' to be used as the adjacent
            %    bulk phase.
            % :parameter area:
            %    The area of the surface in m^2. Defaults to 1.0 m^2 if not
            %    specified.
            % :return:
            %    Instance of class 'ReactorSurface'.

            checklib;

            s.surfID = callct('reactorsurface_new', 0);
            s.reactor = -1;
            if s.surfID < 0
                error(geterr);
            end

            if nargin >= 1
                s.setKinetics(kleft);
            end

            if nargin >= 2
                if isa(reactor, 'Reactor')
                    s.install(reactor);
                else
                    warning('Reactor was not installed due to incorrect type');
                end
            end

            if nargin >= 3
                if isnumeric(area)
                    s.area = area;
                else
                    warning('Area was not a number and was not set');
                end
            end

        end

        %% Utility methods

        function clear(s)
            % Clear the ReactorSurface object from the memory.

            callct('reactorsurface_del', s.surfID);
        end

        function install(s, r)
            % Install a ReactorSurface in a Reactor.

            s.reactor = r;
            callct('reactorsurface_install', s.surfID, r.id);
        end

        function addSensitivityReaction(s, r)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The surface must be
            % installed on a reactor and part of a network first.
            %
            % :parameter m:
            %    Index number of reaction.

            callct('reactorsurface_addSensitivityReaction', s.surfID, r);
        end

        %% ReactorSurface get methods

        function a = get.area(s)
            % Get the areaof the reactor surface in m^2.

            a = callct('reactorsurface_area', s.surfID);
        end

        %% ReactorSurface set methods

        function set.area(s, a)
            % Set the area of a reactor surface

            callct('reactorsurface_setArea', s.surfID, a);
        end

        function setKinetics(s, kin)
            % Setthe surface reaction mechanism on a reactor surface.
            %
            % :parameter kin:
            %    Instance of class 'Kinetics' (or another object derived
            %    from kin) to be used as the kinetic mechanism for this
            %    surface. Typically an instance of class 'Interface'.

            ikin = 0;
            if isa(kin, 'Kinetics')
                ikin = kin.kinID;
            end

            callct('reactorsurface_setkinetics', s.surfID, ikin);
        end
    end
end

