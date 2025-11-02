classdef IdealGasConstPressureReactor < ReactorBase
    % Create a constant pressure reactor with an ideal gas. ::
    %
    %     >> r = IdealGasConstPressureReactor(phase, name, clone)
    %
    % An :mat:class:`IdealGasConstPressureReactor` is an instance of
    % :mat:class:`ReactorBase` where the pressure is held constant.
    % The volume is not a state variable, but instead takes on
    % whatever value is consistent with holding the pressure constant.
    % Additionally, its governing equations are specialized for the
    % ideal gas equation of state (and do not work correctly with other
    % thermodynamic models). Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = IdealGasConstPressureReactor(gas) % a reactor containing a gas
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
    %     Instance of class :mat:class:`IdealGasConstPressureReactor`.

    methods

        function r = IdealGasConstPressureReactor(phase, name, clone)
            % Constructor
            arguments
                phase (1,1) Solution
                name (1,1) string = "(none)"
                clone (1,1) logical = true
            end

            id = ctFunc('mReactor_new', 'IdealGasConstPressureReactor', ...
                        phase.solnID, clone, name);
            r@ReactorBase(id);
        end

    end
end
