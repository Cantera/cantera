classdef ReactorSurface < ReactorBase
    % ReactorSurface Class ::
    %
    %     >> s = ReactorSurface(surf, reactor, name, clone)
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
    % :param reactors:
    %    An instance of or a cell array of instances of class :mat:class:`ReactorBase`.
    % :param name:
    %    Reactor surface name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `surf` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).
    % :return:
    %    Instance of class :mat:class:`ReactorSurface`.

    methods
        %% ReactorSurface Class Constructor

        function s = ReactorSurface(surf, reactors, name, clone)
            % Create a :mat:class:`ReactorSurface` object.

            arguments
                surf {mustBeA(surf, 'Interface')}
                reactors
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            if isa(reactors, 'ReactorBase')
                reactors = {reactors};
            end

            ctIsLoaded;
            reactorIDs = cellfun(@(r) r.id, reactors);
            id = ctFunc('reactor_newSurface', surf.solnID, reactorIDs, clone, name);
            s@ReactorBase(id);
        end

    end

end
