classdef ReactorSurface < handle
    % ReactorSurface Class ::
    %
    %     >> s = ReactorSurface(surf, reactor, area)
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
    % :param surf:
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

        function s = ReactorSurface(surf, reactor, area)
            % Create a :mat:class:`ReactorSurface` object.

            ctIsLoaded;

            s.surfID = ctFunc('reactorsurface_new', 0);
            s.reactor = -1;

            if nargin >= 1
                ikin = 0;

                if isa(surf, 'Kinetics')
                    ikin = surf.kinID;
                end

                ctFunc('reactorsurface_setkinetics', s.surfID, ikin);
            end

            if nargin >= 2

                if isa(reactor, 'Reactor')
                    s.reactor = reactor;
                    ctFunc('reactorsurface_install', s.surfID, reactor.id);
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
            % Delete the :mat:class:`ReactorSurface` object.

            ctFunc('reactorsurface_del', s.surfID);
        end

        %% ReactorSurface Utility Methods

        function addSensitivityReaction(s, m)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The surface must be
            % installed on a reactor and part of a network first.
            %
            % s.addSensitivityReaction(m)
            %
            % :param m:
            %    Index number of reaction.

            ctFunc('reactorsurface_addSensitivityReaction', s.surfID, m);
        end

        %% ReactorSurface Get Methods

        function a = get.area(s)
            a = ctFunc('reactorsurface_area', s.surfID);
        end

        %% ReactorSurface Set Methods

        function set.area(s, a)
            ctFunc('reactorsurface_setArea', s.surfID, a);
        end

    end

end
