classdef ReactorSurface < handle
    % ReactorSurface Class ::
    %
    %     >> s = ReactorSurface(surf, reactor, name)
    %
    % A surface on which heterogeneous reactions take place. The
    % mechanism object (typically an instance of :mat:class:`Interface`)
    % must be constructed so that it is properly linked to the
    % object representing the fluid in the reactor. The surface
    % temperature on each side is taken to be equal to the
    % temperature of the reactor.
    %
    % :param surf:
    %    Surface reaction mechanisms for the left-facing surface.
    %    This must bean instance of class :mat:class:`Kinetics`, or of a class
    %    derived from Kinetics, such as :mat:class:`Interface`.
    % :param reactor:
    %    Instance of class 'Reactor' to be used as the adjacent bulk phase.
    % :param name:
    %    Reactor surface name (optional; default is ``(none)``).
    % :return:
    %    Instance of class :mat:class:`ReactorSurface`.

    properties (SetAccess = immutable)
        surfID
    end

    properties (SetAccess = public)
        name  % Name of reactor surface.

        area % Area of the reactor surface in m^2.
    end

    methods
        %% ReactorSurface Class Constructor

        function s = ReactorSurface(surf, reactor, name)
            % Create a :mat:class:`ReactorSurface` object.

            ctIsLoaded;

            if ~isa(surf, 'Kinetics') || ~isa(reactor, 'Reactor')
                error('Invalid parameters.')
            end
            if nargin < 3
                name = '(none)';
            end

            s.surfID = ctFunc('reactorsurface_new', name);
            ctFunc('reactorsurface_setkinetics', s.surfID, surf.kinID);
            ctFunc('reactorsurface_install', s.surfID, reactor.id);
        end

        %% ReactorSurface Class Destructor

        function delete(s)
            % Delete the :mat:class:`ReactorSurface` object.

            if isempty(s.surfID)
                return
            end
            ctFunc('reactorsurface_del', s.surfID);
        end

        %% ReactorSurface Utility Methods

        function addSensitivityReaction(s, m)
            % Specifies that the sensitivity of the state variables with
            % respect to reaction m should be computed. The surface must be
            % installed on a reactor and part of a network first. ::
            %
            %     >> s.addSensitivityReaction(m)
            %
            % :param m:
            %    Index number of reaction.

            ctFunc('reactorsurface_addSensitivityReaction', s.surfID, m);
        end

        %% ReactorSurface Get Methods

        function name = get.name(s)
            name = ctString('reactorsurface_name', s.surfID);
        end

        function a = get.area(s)
            a = ctFunc('reactorsurface_area', s.surfID);
        end

        %% ReactorSurface Set Methods

        function set.name(s, name)
            ctFunc('reactorsurface_setName', s.surfID, name);
        end

        function set.area(s, a)
            ctFunc('reactorsurface_setArea', s.surfID, a);
        end

    end

end
