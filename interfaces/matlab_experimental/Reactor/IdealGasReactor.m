classdef IdealGasReactor < ReactorBase
    % Create a reactor with an ideal gas. ::
    %
    %     >> r = IdealGasReactor(phase, name, clone)
    %
    % An :mat:class:`IdealGasReactor` is an instance of :mat:class:`ReactorBase` where
    % the governing equations are specialized for the ideal gas equation of state
    % (and do not work correctly with other thermodynamic models). Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = IdealGasReactor(gas)    % a reactor containing a gas
    %
    % See also: :mat:class:`ReactorBase`
    %
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `phase` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).
    % :return:
    %     Instance of class :mat:class:`IdealGasReactor`.

    methods

        function r = IdealGasReactor(phase, name, clone)
            % Constructor

            arguments
                phase {mustBeA(phase, 'Solution')}
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            ctIsLoaded;
            id = ctFunc('reactor_new', 'IdealGasReactor', phase.solnID, clone, name);
            r@ReactorBase(id);
        end

    end

end
