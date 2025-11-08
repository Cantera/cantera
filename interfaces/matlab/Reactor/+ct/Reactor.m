classdef Reactor < ct.ReactorBase
    % Reactor Class ::
    %
    %     >> r = ct.Reactor(phase, name, clone)
    %
    % A :mat:class:`ct.Reactor` object simulates a perfectly-stirred reactor.
    % It has a time-dependent state, and may be coupled to other
    % reactors through flow lines or through walls that may expand
    % or contract and/or conduct heat.
    %
    % :param phase:
    %    Instance of :mat:class:`ct.Solution` representing the contents of
    %    the reactor.
    % :param name:
    %    Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `phase` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).

    methods
        %% Reactor Class Constructor

        function r = Reactor(phase, name, clone)

            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            id = ct.impl.call('mReactor_new', 'Reactor', phase.solnID, clone, name);
            r@ct.ReactorBase(id);
        end

    end

end
