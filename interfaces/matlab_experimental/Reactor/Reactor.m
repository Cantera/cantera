classdef Reactor < ReactorBase
    % Reactor Class ::
    %
    %     >> r = Reactor(content, name, clone)
    %
    % A :mat:class:`Reactor` object simulates a perfectly-stirred reactor.
    % It has a time-dependent state, and may be coupled to other
    % reactors through flow lines or through walls that may expand
    % or contract and/or conduct heat.
    %
    % :param phase:
    %    Instance of :mat:class:`Solution` representing the contents of
    %    the reactor.
    % :param name:
    %    Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `content` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    % :return:
    %    Instance of :mat:class:`Reactor`.

    methods
        %% Reactor Class Constructor

        function r = Reactor(phase, name, clone)

            arguments
                phase {mustBeA(phase, 'Solution')}
                name (1,1) string = "(none)"
                clone (1,1) logical = false
            end

            ctIsLoaded;
            id = ctFunc('reactor_new', 'Reactor', phase.solnID, clone, name);
            r@ReactorBase(id, phase);
        end

    end

end
