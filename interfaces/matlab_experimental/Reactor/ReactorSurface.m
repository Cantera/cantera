classdef ReactorSurface < handle
    % ReactorSurface Class
    %
    % s = ReactorSurface(kleft, reactor, area)
    %
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
    % :param kleft:
    %    Surface reaction mechanisms for the left-facing surface.
    %    This must bean instance of class 'Kinetics', or of a class
    %    derived from 'Kinetics', such as 'Interface'.
    % :param reactor:
    %    Instance of class 'Reactor' to be used as the adjacent
    %    bulk phase.
    % :param area:
    %    The area of the surface in m^2. Defaults to 1.0 m^2 if not
    %    specified.
    % :return:
    %    Instance of class 'ReactorSurface'.
    %

    properties (SetAccess = immutable)
        surfID
    end

    properties (SetAccess = protected)
        reactor
    end

    properties (SetAccess = public)
        area % Area of the reactor surface in m^2.
    end

    methods
        %% ReactorSurface Class Constructor

        function s = ReactorSurface(kleft, reactor, area)
            checklib;

            s.surfID = callct('reactorsurface_new', 0);
            s.reactor = -1;

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

        %% ReactorSurface Class Destructor

        function delete(s)
            % Delete the ReactorSurface object from the memory.

            callct('reactorsurface_del', s.surfID);
        end

        %% ReactorSurface Utility Methods

        function install(s, r)
            % Install a ReactorSurface in a Reactor.
            %
            % s.install(r)
            %
            % :param r:
            %    Instance of class 'Reactor'.
            %

            s.reactor = r;
            callct('reactorsurface_install', s.surfID, r.id);
        end

        function addSensitivityReaction(s, m)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The surface must be
            % installed on a reactor and part of a network first.
            %
            % s.addSensitivityReaction(m)
            %
            % :param m:
            %    Index number of reaction.

            callct('reactorsurface_addSensitivityReaction', s.surfID, m);
        end

        %% ReactorSurface Get Methods

        function a = get.area(s)
            a = callct('reactorsurface_area', s.surfID);
        end

        %% ReactorSurface Set Methods

        function set.area(s, a)
            callct('reactorsurface_setArea', s.surfID, a);
        end

        function setKinetics(s, kin)
            % Setthe surface reaction mechanism on a reactor surface.
            %
            % s.setKinetics(kin)
            %
            % :param kin:
            %    Instance of class 'Kinetics' (or another object derived
            %    from kin) to be used as the kinetic mechanism for this
            %    surface. Typically an instance of class 'Interface'.
            %

            ikin = 0;

            if isa(kin, 'Kinetics')
                ikin = kin.kinID;
            end

            callct('reactorsurface_setkinetics', s.surfID, ikin);
        end

    end

end
