classdef Reservoir < ReactorBase
    % Create a :mat:class:`Reservoir` object. ::
    %
    %     >> r = Reservoir(phase, name)
    %
    % A :mat:class:`Reservoir` is an instance of class :mat:class:`Reactor`
    % configured so that its intensive state is constant in time. A reservoir
    % may be thought of as infinite in extent, perfectly mixed,
    % and non-reacting, so that fluid may be extracted or added without
    % changing the composition or thermodynamic state. Note that even
    % if the reaction mechanism associated with the fluid in the
    % reactor defines reactions, they are disabled within
    % reservoirs. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = Reservoir(gas)    % a reservoir containing a gas
    %
    % See also: :mat:class:`Reactor`
    %
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reservoir name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `content` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    % :return:
    %     Instance of class :mat:class:`Reactor`.

    methods

        function r = Reservoir(phase, name, clone)
            % Constructor

            arguments
                phase {mustBeA(phase, 'Solution')}
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            ctIsLoaded;
            id = ctFunc('reactor_new', 'Reservoir', phase.solnID, clone, name);
            r@ReactorBase(id, phase);
        end

    end
end
