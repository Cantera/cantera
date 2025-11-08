classdef ConstPressureReactor < ct.ReactorBase
    % Create a constant pressure reactor object. ::
    %
    %     >> r = ct.ConstPressureReactor(phase, name, clone)
    %
    % A :mat:class:`ct.ConstPressureReactor` is an instance of class
    % :mat:class:`ct.ReactorBase` where the pressure is held constant. The volume
    % is not a state variable, but instead takes on whatever value is
    % consistent with holding the pressure constant. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = ct.ConstPressureReactor(phase)    % a reactor containing contents
    %
    % See also: :mat:class:`ct.ReactorBase`
    %
    % :param phase:
    %     Cantera :mat:class:`ct.Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :param clone:
    %    Determines whether to clone `phase` so that the internal state of
    %    this reactor is independent of the original Solution object and
    %    any Solution objects used by other reactors in the network.
    %    (optional; default is true).

    methods

        function obj = ConstPressureReactor(phase, name, clone)
            arguments
                phase (1,1) ct.Solution
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            id = ct.impl.call('mReactor_new', 'ConstPressureReactor', ...
                              phase.solnID, clone, name);
            obj@ct.ReactorBase(id);
        end

    end
end
