classdef ReactorSurface < Reactor
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

    properties (SetAccess = public)
        area % Area of the reactor surface in m^2.
    end

    methods
        %% ReactorSurface Class Constructor

        function s = ReactorSurface(surf, reactor, name)
            % Create a :mat:class:`ReactorSurface` object.

            ctIsLoaded;

            if ~isa(surf, 'Solution') || ~isa(reactor, 'Reactor')
                error('Invalid parameters.')
            end
            if nargin < 3
                name = '(none)';
            end

            s@Reactor(surf, 'ReactorSurface', name);
            ctFunc('reactor_addSurface', reactor.id, s.id);
        end

        %% ReactorSurface Get Methods

        function a = get.area(s)
            a = ctFunc('reactor_area', s.id);
        end

        %% ReactorSurface Set Methods

        function set.area(s, a)
            ctFunc('reactor_setArea', s.id, a);
        end

    end

end
