classdef ConstPressureReactor < Reactor
    % Create a constant pressure reactor object. ::
    %
    %     >> r = ConstPressureReactor(phase, name)
    %
    % A :mat:class:`ConstPressureReactor` is an instance of class
    % :mat:class:`Reactor` where the pressure is held constant. The volume
    % is not a state variable, but instead takes on whatever value is
    % consistent with holding the pressure constant. Examples:
    %
    % .. code-block:: matlab
    %
    %     r2 = ConstPressureReactor(phase)    % a reactor containing contents
    %
    % See also: :mat:class:`Reactor`
    %
    % :param phase:
    %     Cantera :mat:class:`Solution` to be set as the contents of the reactor.
    % :param name:
    %     Reactor name (optional; default is ``(none)``).
    % :return:
    %     Instance of class :mat:class:`ConstPressureReactor`.

    methods

        function r = ConstPressureReactor(phase, name)
            % Constructor

            if nargin < 2
                name = '(none)';
            end

            r@Reactor(phase, 'ConstPressureReactor', name);
        end

    end
end
